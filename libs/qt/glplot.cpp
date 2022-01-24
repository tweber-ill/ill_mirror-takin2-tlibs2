/**
 * tlibs2 -- GL plotter
 * @author Tobias Weber <tweber@ill.fr>
 * @date 2017-2021
 * @note The present version was forked on 8-Nov-2018 from my privately developed "magtools" project (https://github.com/t-weber/magtools).
 * @license GPLv3, see 'LICENSE' file
 *
 * References:
 *   - http://doc.qt.io/qt-5/qopenglwidget.html#details
 *   - http://code.qt.io/cgit/qt/qtbase.git/tree/examples/opengl/threadedqopenglwidget
 *
 * ----------------------------------------------------------------------------
 * tlibs
 * Copyright (C) 2017-2021  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * Copyright (C) 2015-2017  Tobias WEBER (Technische Universitaet Muenchen
 *                          (TUM), Garching, Germany).
 * magtools
 * Copyright (C) 2017-2018  Tobias WEBER (privately developed).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * ----------------------------------------------------------------------------
 */

#include "glplot.h"

#include <QtGui/QPainter>
#include <QtGui/QGuiApplication>
#include <QtCore/QtGlobal>

#include <iostream>
#include <boost/scope_exit.hpp>
#include <boost/preprocessor/stringize.hpp>
#include <boost/algorithm/string/replace.hpp>
namespace algo = boost::algorithm;


namespace tl2 {
// ----------------------------------------------------------------------------
// GL plot renderer
// ----------------------------------------------------------------------------

GlPlotRenderer::GlPlotRenderer(GlPlot *pPlot) : m_pPlot{pPlot}
{
	if constexpr(m_usetimer)
	{
		connect(&m_timer, &QTimer::timeout,
			this, static_cast<void (GlPlotRenderer::*)()>(&GlPlotRenderer::tick));
		m_timer.start(std::chrono::milliseconds(1000 / 60));
	}

	UpdateCam();
}


GlPlotRenderer::~GlPlotRenderer()
{
	if constexpr(m_usetimer)
		m_timer.stop();

	// get context
	if constexpr(m_isthreaded)
	{
		m_pPlot->context()->moveToThread(qGuiApp->thread());
	}

	m_pPlot->makeCurrent();
	BOOST_SCOPE_EXIT(m_pPlot) { m_pPlot->doneCurrent(); } BOOST_SCOPE_EXIT_END

	// delete gl objects within current gl context
	m_pShaders.reset();

	qgl_funcs* pGl = GetGlFunctions();
	for(auto &obj : m_objs)
		delete_render_object(obj);

	m_objs.clear();
	LOGGLERR(pGl)
}


void GlPlotRenderer::startedThread() { }
void GlPlotRenderer::stoppedThread() { }


QPointF GlPlotRenderer::GlToScreenCoords(const t_vec_gl& vec4, bool *pVisible) const
{
	t_vec_gl pt = m_cam.ToScreenCoords(vec4, pVisible);
	return QPointF(pt[0], pt[1]);
}


GlPlotObj GlPlotRenderer::CreateTriangleObject(const std::vector<t_vec3_gl>& verts,
	const std::vector<t_vec3_gl>& triagverts, const std::vector<t_vec3_gl>& norms,
	const t_vec_gl& colour, bool bUseVertsAsNorm)
{
	GlPlotObj obj;
	create_triangle_object(m_pPlot, obj, verts, triagverts, norms, {}, colour,
		bUseVertsAsNorm, m_attrVertex, m_attrVertexNorm, m_attrVertexCol, -1);
	return obj;
}


GlPlotObj GlPlotRenderer::CreateLineObject(const std::vector<t_vec3_gl>& verts, const t_vec_gl& colour)
{
	GlPlotObj obj;
	create_line_object(m_pPlot, obj, verts, colour, m_attrVertex, m_attrVertexCol);
	return obj;
}


void GlPlotRenderer::SetObjectMatrix(std::size_t idx, const t_mat_gl& mat)
{
	if(idx >= m_objs.size()) return;
	m_objs[idx].m_mat = mat;
}


void GlPlotRenderer::SetObjectCol(std::size_t idx, t_real_gl r, t_real_gl g, t_real_gl b, t_real_gl a)
{
	if(idx >= m_objs.size()) return;
	m_objs[idx].m_colour = tl2::create<t_vec_gl>({r,g,b,a});
}


void GlPlotRenderer::SetObjectLabel(std::size_t idx, const std::string& label)
{
	if(idx >= m_objs.size()) return;
	m_objs[idx].m_label = label;
}

const std::string& GlPlotRenderer::GetObjectLabel(std::size_t idx) const
{
	static const std::string empty{};
	if(idx >= m_objs.size()) return empty;

	return m_objs[idx].m_label;
}


void GlPlotRenderer::SetObjectDataString(std::size_t idx, const std::string& data)
{
	if(idx >= m_objs.size()) return;
	m_objs[idx].m_datastr = data;
}

const std::string& GlPlotRenderer::GetObjectDataString(std::size_t idx) const
{
	static const std::string empty{};
	if(idx >= m_objs.size()) return empty;

	return m_objs[idx].m_datastr;
}


void GlPlotRenderer::SetObjectVisible(std::size_t idx, bool visible)
{
	if(idx >= m_objs.size()) return;
	m_objs[idx].m_visible = visible;
}


void GlPlotRenderer::SetObjectHighlight(std::size_t idx, bool highlight)
{
	if(idx >= m_objs.size()) return;
	m_objs[idx].m_highlighted = highlight;
}


bool GlPlotRenderer::GetObjectHighlight(std::size_t idx) const
{
	if(idx >= m_objs.size()) return 0;
	return m_objs[idx].m_highlighted;
}


void GlPlotRenderer::RemoveObject(std::size_t idx)
{
	m_objs[idx].m_valid = false;

	m_objs[idx].m_vertex_buffer.reset();
	m_objs[idx].m_normals_buffer.reset();
	m_objs[idx].m_colour_buffer.reset();

	m_objs[idx].m_vertices.clear();
	m_objs[idx].m_triangles.clear();

	// TODO: remove if object has no follow-up indices
}


std::size_t GlPlotRenderer::AddLinkedObject(std::size_t linkTo,
	t_real_gl x, t_real_gl y, t_real_gl z,
	t_real_gl r, t_real_gl g, t_real_gl b, t_real_gl a)
{
	GlPlotObj obj;
	obj.linkedObj = linkTo;
	obj.m_mat = tl2::hom_translation<t_mat_gl>(x, y, z);
	obj.m_colour = tl2::create<t_vec_gl>({r, g, b, a});

	QMutexLocker _locker{&m_mutexObj};
	m_objs.emplace_back(std::move(obj));

	return m_objs.size()-1;		// object handle
}


std::size_t GlPlotRenderer::AddSphere(t_real_gl rad, t_real_gl x, t_real_gl y, t_real_gl z,
	t_real_gl r, t_real_gl g, t_real_gl b, t_real_gl a)
{
	constexpr int numsubdivs = 2;

	auto solid = tl2::create_icosahedron<t_vec3_gl>(1);
	auto [triagverts, norms, uvs] = tl2::spherify<t_vec3_gl>(
		tl2::subdivide_triangles<t_vec3_gl>(
			tl2::create_triangles<t_vec3_gl>(solid), numsubdivs), rad);
	auto [boundingSpherePos, boundingSphereRad] = tl2::bounding_sphere<t_vec3_gl>(triagverts);

	QMutexLocker _locker{&m_mutexObj};

	auto obj = CreateTriangleObject(std::get<0>(solid), triagverts, norms, tl2::create<t_vec_gl>({r,g,b,a}), true);
	obj.m_mat = tl2::hom_translation<t_mat_gl>(x, y, z);
	obj.m_boundingSpherePos = std::move(boundingSpherePos);
	obj.m_boundingSphereRad = boundingSphereRad;
	//obj.m_boundingSphereRad = rad;
	m_objs.emplace_back(std::move(obj));

	return m_objs.size()-1;		// object handle
}


std::size_t GlPlotRenderer::AddCylinder(t_real_gl rad, t_real_gl h,
	t_real_gl x, t_real_gl y, t_real_gl z,
	t_real_gl r, t_real_gl g, t_real_gl b, t_real_gl a)
{
	auto solid = tl2::create_cylinder<t_vec3_gl>(rad, h, true);
	auto [triagverts, norms, uvs] = tl2::create_triangles<t_vec3_gl>(solid);
	auto [boundingSpherePos, boundingSphereRad] = tl2::bounding_sphere<t_vec3_gl>(triagverts);

	QMutexLocker _locker{&m_mutexObj};

	auto obj = CreateTriangleObject(std::get<0>(solid), triagverts, norms, tl2::create<t_vec_gl>({r,g,b,a}), false);
	obj.m_mat = tl2::hom_translation<t_mat_gl>(x, y, z);
	obj.m_boundingSpherePos = std::move(boundingSpherePos);
	obj.m_boundingSphereRad = boundingSphereRad;
	m_objs.emplace_back(std::move(obj));

	return m_objs.size()-1;		// object handle
}


std::size_t GlPlotRenderer::AddCone(t_real_gl rad, t_real_gl h,
	t_real_gl x, t_real_gl y, t_real_gl z,
	t_real_gl r, t_real_gl g, t_real_gl b, t_real_gl a)
{
	auto solid = tl2::create_cone<t_vec3_gl>(rad, h);
	auto [triagverts, norms, uvs] = tl2::create_triangles<t_vec3_gl>(solid);
	auto [boundingSpherePos, boundingSphereRad] = tl2::bounding_sphere<t_vec3_gl>(triagverts);

	QMutexLocker _locker{&m_mutexObj};

	auto obj = CreateTriangleObject(std::get<0>(solid), triagverts, norms, tl2::create<t_vec_gl>({r,g,b,a}), false);
	obj.m_mat = tl2::hom_translation<t_mat_gl>(x, y, z);
	obj.m_boundingSpherePos = std::move(boundingSpherePos);
	obj.m_boundingSphereRad = boundingSphereRad;
	m_objs.emplace_back(std::move(obj));

	return m_objs.size()-1;		// object handle
}


std::size_t GlPlotRenderer::AddArrow(t_real_gl rad, t_real_gl h,
	t_real_gl x, t_real_gl y, t_real_gl z,
	t_real_gl r, t_real_gl g, t_real_gl b, t_real_gl a)
{
	auto solid = tl2::create_cylinder<t_vec3_gl>(rad, h, 2, 32, rad, rad*1.5);
	auto [triagverts, norms, uvs] = tl2::create_triangles<t_vec3_gl>(solid);
	auto [boundingSpherePos, boundingSphereRad] = tl2::bounding_sphere<t_vec3_gl>(triagverts);

	QMutexLocker _locker{&m_mutexObj};

	auto obj = CreateTriangleObject(std::get<0>(solid), triagverts, norms, tl2::create<t_vec_gl>({r,g,b,a}), false);
	obj.m_mat = tl2::get_arrow_matrix<t_vec_gl, t_mat_gl, t_real_gl>(
		tl2::create<t_vec_gl>({1,0,0}), 1., tl2::create<t_vec_gl>({x,y,z}), tl2::create<t_vec_gl>({0,0,1}));
	obj.m_boundingSpherePos = std::move(boundingSpherePos);
	obj.m_boundingSphereRad = boundingSphereRad;
	obj.m_labelPos = tl2::create<t_vec3_gl>({0., 0., 0.75});
	m_objs.emplace_back(std::move(obj));

	return m_objs.size()-1;		// object handle
}


std::size_t GlPlotRenderer::AddTriangleObject(const std::vector<t_vec3_gl>& triag_verts,
	const std::vector<t_vec3_gl>& triag_norms,
	t_real_gl r, t_real_gl g, t_real_gl b, t_real_gl a)
{
	auto [boundingSpherePos, boundingSphereRad] = tl2::bounding_sphere<t_vec3_gl>(triag_verts);

	QMutexLocker _locker{&m_mutexObj};

	auto obj = CreateTriangleObject(triag_verts, triag_verts, triag_norms, tl2::create<t_vec_gl>({r,g,b,a}), false);
	obj.m_mat = tl2::hom_translation<t_mat_gl, t_real_gl>(0., 0., 0.);
	obj.m_boundingSpherePos = std::move(boundingSpherePos);
	obj.m_boundingSphereRad = boundingSphereRad;
	obj.m_labelPos = tl2::create<t_vec3_gl>({0., 0., 0.75});
	m_objs.emplace_back(std::move(obj));

	return m_objs.size()-1;		// object handle
}


std::size_t GlPlotRenderer::AddCoordinateCross(t_real_gl min, t_real_gl max)
{
	auto col = tl2::create<t_vec_gl>({0,0,0,1});
	auto verts = std::vector<t_vec3_gl>
	{{
		tl2::create<t_vec3_gl>({min,0,0}), tl2::create<t_vec3_gl>({max,0,0}),
		tl2::create<t_vec3_gl>({0,min,0}), tl2::create<t_vec3_gl>({0,max,0}),
		tl2::create<t_vec3_gl>({0,0,min}), tl2::create<t_vec3_gl>({0,0,max}),
	}};

	QMutexLocker _locker{&m_mutexObj};

	auto obj = CreateLineObject(verts, col);
	obj.m_invariant = true;
	m_objs.emplace_back(std::move(obj));

	return m_objs.size()-1;		// object handle
}



void GlPlotRenderer::initialiseGL()
{
	// --------------------------------------------------------------------
	// shaders
	// --------------------------------------------------------------------
	std::string strFragShader = R"RAW(#version ${GLSL_VERSION}

in vec4 fragcol;
out vec4 outcol;


void main()
{
	outcol = fragcol;
})RAW";
	// --------------------------------------------------------------------


	// --------------------------------------------------------------------
	std::string strVertexShader = R"RAW(#version ${GLSL_VERSION}

in vec4 vertex;
in vec4 normal;
in vec4 vertexcol;
out vec4 fragcol;


const float pi = ${PI};


// ----------------------------------------------------------------------------
// transformations
// ----------------------------------------------------------------------------
uniform mat4 proj = mat4(1.);
uniform mat4 cam = mat4(1.);
uniform mat4 cam_inv = mat4(1.);
uniform mat4 obj = mat4(1.);
uniform mat4 trafoA = mat4(1.);
uniform mat4 trafoB = mat4(1.);		// B = 2 pi / A

uniform int coordsys = 0;			// 0: crystal system, 1: lab system
// ----------------------------------------------------------------------------


// ----------------------------------------------------------------------------
// lighting
// ----------------------------------------------------------------------------
uniform vec4 constcol = vec4(1, 1, 1, 1);
uniform vec3 lightpos[] = vec3[]( vec3(5, 5, 5), vec3(0, 0, 0), vec3(0, 0, 0), vec3(0, 0, 0) );
uniform int activelights = 1;	// how many lights to use?

float g_diffuse = 1.;
float g_specular = 0.25;
float g_shininess = 1.;
float g_ambient = 0.2;
// ----------------------------------------------------------------------------


/**
 * reflect a vector on a surface with normal n => subtract the projection vector twice: 1 - 2*|n><n|
 */
mat3 reflect(vec3 n)
{
	mat3 refl = mat3(1.) - 2.*outerProduct(n, n);

	// have both vectors point away from the surface
	return -refl;
}


/**
 * position of the camera
 */
vec3 get_campos()
{
	vec4 trans = -vec4(cam[3].xyz, 0);
	return (cam_inv*trans).xyz;
}


/**
 * phong lighting model
 * @see: https://en.wikipedia.org/wiki/Phong_reflection_model
 */
float lighting(vec4 objVert, vec4 objNorm)
{
	float I_diff = 0.;
	float I_spec = 0.;


	vec3 dirToCam;
	// only used for specular lighting
	if(g_specular > 0.) dirToCam = normalize(get_campos() - objVert.xyz);


	// iterate (active) light sources
	for(int lightidx=0; lightidx<min(lightpos.length(), activelights); ++lightidx)
	{
		// diffuse lighting
		vec3 dirLight = normalize(lightpos[lightidx]-objVert.xyz);

		if(g_diffuse > 0.)
		{
			float I_diff_inc = g_diffuse * dot(objNorm.xyz, dirLight);
			if(I_diff_inc < 0.) I_diff_inc = 0.;
			I_diff += I_diff_inc;
		}


		// specular lighting
		if(g_specular > 0.)
		{
			if(dot(dirToCam, objNorm.xyz) > 0.)
			{
				vec3 dirLightRefl = reflect(objNorm.xyz) * dirLight;

				float val = dot(dirToCam, dirLightRefl);
				if(val > 0.)
				{
					float I_spec_inc = g_specular * pow(val, g_shininess);
					if(I_spec_inc < 0.) I_spec_inc = 0.;
					I_spec += I_spec_inc;
				}
			}
		}
	}


	// ambient lighting
	float I_amb = g_ambient;


	// total intensity
	return I_diff + I_spec + I_amb;
}


void main()
{
	mat4 coordTrafo = mat4(1.);
	mat4 coordTrafo_inv = mat4(1.);

	if(coordsys == 1)
	{
		coordTrafo = trafoA;
		coordTrafo_inv = trafoB / (2.*pi);
		coordTrafo_inv[3][3] = 1.;
	}

	// coordTrafo_inv is needed so not to distort the object
	vec4 objPos = coordTrafo * obj * coordTrafo_inv * vertex;
	vec4 objNorm = normalize(coordTrafo * obj * coordTrafo_inv * normal);
	gl_Position = proj * cam * objPos;

	float I = lighting(objPos, objNorm);
	fragcol.rgb = vertexcol.rgb * I;
	fragcol *= constcol;
})RAW";
// --------------------------------------------------------------------


	// set glsl version and constants
	const std::string strGlsl = std::to_string(_GLSL_MAJ_VER*100 + _GLSL_MIN_VER*10);
	std::string strPi = std::to_string(tl2::pi<t_real_gl>);	// locale-dependent !
	algo::replace_all(strPi, std::string(","), std::string("."));	// ensure decimal point

	for(std::string* strSrc : { &strFragShader, &strVertexShader })
	{
		algo::replace_all(*strSrc, std::string("${GLSL_VERSION}"), strGlsl);
		algo::replace_all(*strSrc, std::string("${PI}"), strPi);
	}


	// GL functions
	auto *pGl = GetGlFunctions();
	if(!pGl) return;

	m_strGlVer = (char*)pGl->glGetString(GL_VERSION);
	m_strGlShaderVer = (char*)pGl->glGetString(GL_SHADING_LANGUAGE_VERSION);
	m_strGlVendor = (char*)pGl->glGetString(GL_VENDOR);
	m_strGlRenderer = (char*)pGl->glGetString(GL_RENDERER);
	LOGGLERR(pGl);


	// shaders
	{
		static QMutex shadermutex;
		shadermutex.lock();
		BOOST_SCOPE_EXIT(&shadermutex) { shadermutex.unlock(); } BOOST_SCOPE_EXIT_END

		// shader compiler/linker error handler
		auto shader_err = [this](const char* err) -> void
		{
			std::cerr << err << std::endl;

			std::string strLog = m_pShaders->log().toStdString();
			if(strLog.size())
				std::cerr << "Shader log: " << strLog << std::endl;

			std::exit(-1);
		};

		// compile & link shaders
		m_pShaders = std::make_shared<QOpenGLShaderProgram>(this);

		if(!m_pShaders->addShaderFromSourceCode(QOpenGLShader::Fragment, strFragShader.c_str()))
			shader_err("Cannot compile fragment shader.");
		if(!m_pShaders->addShaderFromSourceCode(QOpenGLShader::Vertex, strVertexShader.c_str()))
			shader_err("Cannot compile vertex shader.");

		if(!m_pShaders->link())
			shader_err("Cannot link shaders.");

		m_uniMatrixCam = m_pShaders->uniformLocation("cam");
		m_uniMatrixCamInv = m_pShaders->uniformLocation("cam_inv");
		m_uniMatrixProj = m_pShaders->uniformLocation("proj");
		m_uniMatrixObj = m_pShaders->uniformLocation("obj");
		m_uniMatrixA = m_pShaders->uniformLocation("trafoA");
		m_uniMatrixB = m_pShaders->uniformLocation("trafoB");
		m_uniCoordSys = m_pShaders->uniformLocation("coordsys");
		m_uniConstCol = m_pShaders->uniformLocation("constcol");
		m_uniLightPos = m_pShaders->uniformLocation("lightpos");
		m_uniNumActiveLights = m_pShaders->uniformLocation("activelights");
		m_attrVertex = m_pShaders->attributeLocation("vertex");
		m_attrVertexNorm = m_pShaders->attributeLocation("normal");
		m_attrVertexCol = m_pShaders->attributeLocation("vertexcol");
	}
	LOGGLERR(pGl);


	// 3d objects
	AddCoordinateCross(-m_CoordMax, m_CoordMax);


	m_bInitialised = true;

	// check threading compatibility
	if constexpr(m_isthreaded)
	{
		if(auto *pContext = ((QOpenGLWidget*)m_pPlot)->context(); pContext && !pContext->supportsThreadedOpenGL())
		{
			m_bPlatformSupported = false;
			std::cerr << "Threaded GL is not supported on this platform." << std::endl;
		}
	}
}


void GlPlotRenderer::SetScreenDims(int w, int h)
{
	m_cam.SetScreenDimensions(w, h);
	m_bWantsResize = true;
}


void GlPlotRenderer::resizeGL()
{
	if(!m_bPlatformSupported || !m_bInitialised) return;

	const auto [w, h] = m_cam.GetScreenDimensions();
	const auto [z_near, z_far] = m_cam.GetDepthRange();

	if(auto *pContext = ((QOpenGLWidget*)m_pPlot)->context(); !pContext)
		return;
	auto *pGl = GetGlFunctions();
	if(!pGl)
		return;

	m_cam.UpdateViewport();
	m_cam.UpdatePerspective();

	pGl->glViewport(0, 0, w, h);
	pGl->glDepthRange(z_near, z_far);

	// bind shaders
	m_pShaders->bind();
	BOOST_SCOPE_EXIT(m_pShaders) { m_pShaders->release(); } BOOST_SCOPE_EXIT_END
	LOGGLERR(pGl);

	// set matrices
	m_pShaders->setUniformValue(m_uniMatrixCam, m_cam.GetTransformation());
	m_pShaders->setUniformValue(m_uniMatrixCamInv, m_cam.GetInverseTransformation());
	m_pShaders->setUniformValue(m_uniMatrixProj, m_cam.GetPerspective());
	LOGGLERR(pGl);

	m_bWantsResize = false;
}


/**
 * set up a (crystal) B matrix
 */
void GlPlotRenderer::SetBTrafo(const t_mat_gl& matB, const t_mat_gl* matA)
{
	m_matB = matB;

	// if A matix is not given, calculate it
	if(matA)
	{
		m_matA = *matA;
	}
	else
	{
		bool ok = true;
		std::tie(m_matA, ok) = tl2::inv(m_matB);
		if(!ok)
		{
			m_matA = tl2::unit<t_mat_gl>();
			std::cerr << "Error: Cannot invert B matrix." << std::endl;
		}
		else
		{
			m_matA *= t_real_gl(2)*tl2::pi<t_real_gl>;
			m_matA(3,3) = 1;
		}
	}

	m_bBTrafoNeedsUpdate = true;
	RequestPlotUpdate();
}


void GlPlotRenderer::SetCoordSys(int iSys)
{
	m_iCoordSys = iSys;
	RequestPlotUpdate();
}


/**
 * update the shader's B matrix
 */
void GlPlotRenderer::UpdateBTrafo()
{
	m_pShaders->setUniformValue(m_uniMatrixA, m_matA);
	m_pShaders->setUniformValue(m_uniMatrixB, m_matB);

	m_bBTrafoNeedsUpdate = false;
}


void GlPlotRenderer::UpdateCam()
{
	m_cam.UpdateTransformation();

	m_bPickerNeedsUpdate = true;
	RequestPlotUpdate();
}


/**
 * request a plot update
 */
void GlPlotRenderer::RequestPlotUpdate()
{
#if QT_VERSION >= QT_VERSION_CHECK(5, 10, 0)
	QMetaObject::invokeMethod((QOpenGLWidget*)m_pPlot,
		static_cast<void (QOpenGLWidget::*)()>(&QOpenGLWidget::update),
		Qt::ConnectionType::QueuedConnection);
#else
	QMetaObject::invokeMethod((QOpenGLWidget*)m_pPlot,
		"update",
		Qt::ConnectionType::QueuedConnection);
#endif
}


void GlPlotRenderer::SetLight(std::size_t idx, const t_vec3_gl& pos)
{
	if(m_lights.size() < idx+1)
		m_lights.resize(idx+1);

	m_lights[idx] = pos;
	m_bLightsNeedUpdate = true;
}


void GlPlotRenderer::UpdateLights()
{
	constexpr int MAX_LIGHTS = 4;	// max. number allowed in shader

	int num_lights = std::min(MAX_LIGHTS, static_cast<int>(m_lights.size()));
	t_real_gl pos[num_lights * 3];

	for(int i=0; i<num_lights; ++i)
	{
		pos[i*3 + 0] = m_lights[i][0];
		pos[i*3 + 1] = m_lights[i][1];
		pos[i*3 + 2] = m_lights[i][2];
	}

	m_pShaders->setUniformValueArray(m_uniLightPos, pos, num_lights, 3);
	m_pShaders->setUniformValue(m_uniNumActiveLights, num_lights);

	m_bLightsNeedUpdate = false;
}


void GlPlotRenderer::EnablePicker(bool b)
{
	m_bPickerEnabled = b;
}


void GlPlotRenderer::UpdatePicker()
{
	if(!m_bInitialised || !m_bPlatformSupported || !m_bPickerEnabled) return;

	// picker ray
	auto [org3, dir3] = m_cam.GetPickerRay(m_posMouse.x(), m_posMouse.y());

	// intersection with unit sphere around origin
	bool hasSphereInters = false;
	t_vec_gl vecClosestSphereInters = tl2::create<t_vec_gl>({0,0,0,0});

	auto intersUnitSphere =
	tl2::intersect_line_sphere<t_vec3_gl, std::vector>(org3, dir3,
		tl2::create<t_vec3_gl>({0,0,0}), t_real_gl(m_pickerSphereRadius));
	for(const auto& result : intersUnitSphere)
	{
		t_vec_gl vecInters4 = tl2::create<t_vec_gl>({result[0], result[1], result[2], 1});

		if(!hasSphereInters)
		{	// first intersection
			vecClosestSphereInters = vecInters4;
			hasSphereInters = true;
		}
		else
		{	// test if next intersection is closer...
			t_vec_gl oldPosTrafo = m_cam.GetTransformation() * vecClosestSphereInters;
			t_vec_gl newPosTrafo = m_cam.GetTransformation() * vecInters4;

			// ... it is closer.
			if(tl2::norm(newPosTrafo) < tl2::norm(oldPosTrafo))
				vecClosestSphereInters = vecInters4;
		}
	}


	// crystal or lab coordinate system?
	const t_mat_gl matUnit = tl2::unit<t_mat_gl>();
	const t_mat_gl *coordTrafo = &matUnit;
	t_mat_gl coordTrafoInv = matUnit;
	if(m_iCoordSys == 1)
	{
		coordTrafo = &m_matA;
		coordTrafoInv = m_matB / (t_real_gl(2)*tl2::pi<t_real_gl>);
		coordTrafoInv(3,3) = 1;
	}


	// intersection with geometry
	bool hasInters = false;
	t_vec_gl vecClosestInters = tl2::create<t_vec_gl>({0,0,0,0});
	std::size_t objInters = 0xffffffff;


	QMutexLocker _locker{&m_mutexObj};

	for(std::size_t curObj=0; curObj<m_objs.size(); ++curObj)
	{
		const auto& obj = m_objs[curObj];
		const GlPlotObj *linkedObj = &obj;
		if(obj.linkedObj)
			linkedObj = &m_objs[*obj.linkedObj];

		if(linkedObj->m_type != GlPlotObjType::TRIANGLES || !obj.m_visible || !obj.m_valid)
			continue;


		const t_mat_gl& matTrafo = (*coordTrafo) * obj.m_mat * coordTrafoInv;

		// scaling factor, TODO: maximum factor for non-uniform scaling
		auto scale = std::cbrt(std::abs(tl2::det(matTrafo)));

		// intersection with bounding sphere?
		auto boundingInters =
			tl2::intersect_line_sphere<t_vec3_gl, std::vector>(org3, dir3,
				matTrafo * linkedObj->m_boundingSpherePos, scale*linkedObj->m_boundingSphereRad);
		if(boundingInters.size() == 0)
			continue;


		// test actual polygons for intersection
		for(std::size_t startidx=0; startidx+2<linkedObj->m_triangles.size(); startidx+=3)
		{
			std::vector<t_vec3_gl> poly{ {
				linkedObj->m_triangles[startidx+0],
				linkedObj->m_triangles[startidx+1],
				linkedObj->m_triangles[startidx+2]
			} };

			/*std::vector<t_vec3_gl> polyuv{ { 
				linkedObj->m_uvs[startidx+0], 
				linkedObj->m_uvs[startidx+1], 
				linkedObj->m_uvs[startidx+2]
			} };*/


			// coordTrafoInv only keeps 3d objects from locally distorting
			auto [vecInters, bInters, lamInters] =
				tl2::intersect_line_poly<t_vec3_gl, t_mat_gl>(org3, dir3, poly, matTrafo);

			if(bInters)
			{
				t_vec_gl vecInters4 = tl2::create<t_vec_gl>({vecInters[0], vecInters[1], vecInters[2], 1});

				if(!hasInters)
				{	// first intersection
					vecClosestInters = vecInters4;
					objInters = curObj;
					hasInters = true;
				}
				else
				{	// test if next intersection is closer...
					t_vec_gl oldPosTrafo = m_cam.GetTransformation() * vecClosestInters;
					t_vec_gl newPosTrafo = m_cam.GetTransformation() * vecInters4;

					if(tl2::norm(newPosTrafo) < tl2::norm(oldPosTrafo))
					{	// ...it is closer
						vecClosestInters = vecInters4;
						objInters = curObj;
					}
				}

				// intersection point in uv coordinates:
				//auto uv = tl2::poly_uv<t_mat_gl, t_vec3_gl>
				//	(poly[0], poly[1], poly[2], polyuv[0], polyuv[1], polyuv[2], vecInters);
			}
		}
	}

	m_bPickerNeedsUpdate = false;
	t_vec3_gl vecClosestInters3 = tl2::create<t_vec3_gl>({vecClosestInters[0], vecClosestInters[1], vecClosestInters[2]});
	t_vec3_gl vecClosestSphereInters3 = tl2::create<t_vec3_gl>({vecClosestSphereInters[0], vecClosestSphereInters[1], vecClosestSphereInters[2]});
	emit PickerIntersection(hasInters ? &vecClosestInters3 : nullptr, objInters,
		hasSphereInters ? &vecClosestSphereInters3 : nullptr);
}


void GlPlotRenderer::mouseMoveEvent(const QPointF& pos)
{
	m_posMouse = pos;

	if(m_bInRotation)
	{
		auto diff = m_posMouse - m_posMouseRotationStart;

		m_cam.Rotate(
			diff.x() / 180. * tl2::pi<t_real_gl>,
			diff.y() / 180. * tl2::pi<t_real_gl>);
		UpdateCam();
	}
	else
	{
		// also automatically done in UpdateCam
		m_bPickerNeedsUpdate = true;
		RequestPlotUpdate();
	}
}


void GlPlotRenderer::zoom(t_real_gl val)
{
	m_cam.Zoom(val/64.);
	UpdateCam();
}


void GlPlotRenderer::ResetZoom()
{
	m_cam.SetZoom(1.);
	UpdateCam();
}


void GlPlotRenderer::BeginRotation()
{
	if(!m_bInRotation)
	{
		m_posMouseRotationStart = m_posMouse;
		m_bInRotation = true;
	}
}


void GlPlotRenderer::EndRotation()
{
	if(m_bInRotation)
	{
		m_cam.SaveRotation();
		m_bInRotation = false;
	}
}


void GlPlotRenderer::tick()
{
	tick(std::chrono::milliseconds(1000 / 60));
}


void GlPlotRenderer::tick(const std::chrono::milliseconds& ms)
{
	// TODO
	UpdateCam();
}


/**
 * pure gl drawing
 */
void GlPlotRenderer::DoPaintGL(qgl_funcs *pGl)
{
	if(!m_bInitialised || !pGl || thread() != QThread::currentThread())
		return;

	// options
	pGl->glCullFace(GL_BACK);
	pGl->glFrontFace(GL_CCW);
	pGl->glEnable(GL_CULL_FACE);

	pGl->glDisable(GL_BLEND);

	pGl->glEnable(GL_MULTISAMPLE);
	pGl->glEnable(GL_LINE_SMOOTH);
	pGl->glEnable(GL_POLYGON_SMOOTH);
	pGl->glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
	pGl->glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);

	// clear
	pGl->glClearColor(1., 1., 1., 1.);
	pGl->glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	pGl->glEnable(GL_DEPTH_TEST);
	LOGGLERR(pGl);


	// bind shaders
	m_pShaders->bind();
	BOOST_SCOPE_EXIT(m_pShaders) { m_pShaders->release(); } BOOST_SCOPE_EXIT_END
	LOGGLERR(pGl);

	if(m_bLightsNeedUpdate) UpdateLights();
	if(m_bBTrafoNeedsUpdate) UpdateBTrafo();

	// set cam matrix
	m_pShaders->setUniformValue(m_uniMatrixCam, m_cam.GetTransformation());
	m_pShaders->setUniformValue(m_uniMatrixCamInv, m_cam.GetInverseTransformation());


	auto colOverride = tl2::create<t_vec_gl>({1,1,1,1});
	auto colHighlight = tl2::create<t_vec_gl>({1,1,1,1});

	// render triangle geometry
	for(const auto& obj : m_objs)
	{
		const GlPlotObj *linkedObj = &obj;
		if(obj.linkedObj)
		{
			// get linked object
			linkedObj = &m_objs[*obj.linkedObj];

			// override constant colour for linked object
			if(obj.m_highlighted)
				m_pShaders->setUniformValue(m_uniConstCol, colHighlight);
			else
				m_pShaders->setUniformValue(m_uniConstCol, obj.m_colour);
		}
		else
		{
			// set override colour to white for non-linked objects
			m_pShaders->setUniformValue(m_uniConstCol, colOverride);
		}

		if(!obj.m_visible || !obj.m_valid) continue;


		m_pShaders->setUniformValue(m_uniMatrixObj, obj.m_mat);

		// set to untransformed coordinate system if the object is invariant
		m_pShaders->setUniformValue(m_uniCoordSys, linkedObj->m_invariant ? 0 : m_iCoordSys.load());


		// main vertex array object
		linkedObj->m_vertex_array->bind();

		pGl->glEnableVertexAttribArray(m_attrVertex);
		if(linkedObj->m_type == GlPlotObjType::TRIANGLES)
			pGl->glEnableVertexAttribArray(m_attrVertexNorm);
		pGl->glEnableVertexAttribArray(m_attrVertexCol);
		BOOST_SCOPE_EXIT(pGl, &m_attrVertex, &m_attrVertexNorm, &m_attrVertexCol)
		{
			pGl->glDisableVertexAttribArray(m_attrVertexCol);
			pGl->glDisableVertexAttribArray(m_attrVertexNorm);
			pGl->glDisableVertexAttribArray(m_attrVertex);
		}
		BOOST_SCOPE_EXIT_END
		LOGGLERR(pGl);


		if(linkedObj->m_type == GlPlotObjType::TRIANGLES)
			pGl->glDrawArrays(GL_TRIANGLES, 0, linkedObj->m_triangles.size());
		else if(linkedObj->m_type == GlPlotObjType::LINES)
			pGl->glDrawArrays(GL_LINES, 0, linkedObj->m_vertices.size());
		else
			std::cerr << "Unknown plot object type." << std::endl;

		LOGGLERR(pGl);
	}

	pGl->glDisable(GL_DEPTH_TEST);
}


/**
 * directly draw on a qpainter
 */
void GlPlotRenderer::DoPaintNonGL(QPainter &painter)
{
	const t_mat_gl matUnit = tl2::unit<t_mat_gl>();

	QFont fontOrig = painter.font();
	QPen penOrig = painter.pen();

	QPen penLabel(Qt::black);
	painter.setPen(penLabel);


	// coordinate labels
	painter.drawText(GlToScreenCoords(tl2::create<t_vec_gl>({0.,0.,0.,1.})), "0");
	for(t_real_gl f=-std::floor(m_CoordMax); f<=std::floor(m_CoordMax); f+=0.5)
	{
		if(tl2::equals<t_real_gl>(f, 0))
			continue;

		std::ostringstream ostrF;
		ostrF << f;
		painter.drawText(GlToScreenCoords(tl2::create<t_vec_gl>({f,0.,0.,1.})), ostrF.str().c_str());
		painter.drawText(GlToScreenCoords(tl2::create<t_vec_gl>({0.,f,0.,1.})), ostrF.str().c_str());
		painter.drawText(GlToScreenCoords(tl2::create<t_vec_gl>({0.,0.,f,1.})), ostrF.str().c_str());
	}

	painter.drawText(GlToScreenCoords(tl2::create<t_vec_gl>({m_CoordMax*t_real_gl(1.2), 0., 0., 1.})), "x");
	painter.drawText(GlToScreenCoords(tl2::create<t_vec_gl>({0., m_CoordMax*t_real_gl(1.2), 0., 1.})), "y");
	painter.drawText(GlToScreenCoords(tl2::create<t_vec_gl>({0., 0., m_CoordMax*t_real_gl(1.2), 1.})), "z");


	// render object labels
	for(const auto& obj : m_objs)
	{
		if(!obj.m_visible || !obj.m_valid) continue;

		if(obj.m_label != "")
		{
			const t_mat_gl *coordTrafo = &matUnit;
			if(m_iCoordSys == 1 && !obj.m_invariant) coordTrafo = &m_matA;

			t_vec3_gl posLabel3d = (*coordTrafo) * obj.m_mat * obj.m_labelPos;
			auto posLabel2d = GlToScreenCoords(tl2::create<t_vec_gl>({posLabel3d[0], posLabel3d[1], posLabel3d[2], 1.}));

			QFont fontLabel = fontOrig;
			QPen penLabel = penOrig;

			fontLabel.setStyleStrategy(QFont::StyleStrategy(/*QFont::OpenGLCompatible |*/ QFont::PreferAntialias | QFont::PreferQuality));
			fontLabel.setWeight(QFont::Medium);
			//penLabel.setColor(QColor(int((1.-obj.m_colour[0])*255.), int((1.-obj.m_colour[1])*255.), int((1.-obj.m_colour[2])*255.), int(obj.m_colour[3]*255.)));
			penLabel.setColor(QColor(0,0,0,255));
			painter.setFont(fontLabel);
			painter.setPen(penLabel);
			painter.drawText(posLabel2d, obj.m_label.c_str());

			fontLabel.setWeight(QFont::Normal);
			penLabel.setColor(QColor(int(obj.m_colour[0]*255.), int(obj.m_colour[1]*255.), int(obj.m_colour[2]*255.), int(obj.m_colour[3]*255.)));
			painter.setFont(fontLabel);
			painter.setPen(penLabel);
			painter.drawText(posLabel2d, obj.m_label.c_str());
		}
	}

	// restore original styles
	painter.setFont(fontOrig);
	painter.setPen(penOrig);
}


void GlPlotRenderer::paintGL()
{
	if(!m_bPlatformSupported) return;
	QMutexLocker _locker{&m_mutexObj};

	if constexpr(!m_isthreaded)
	{
		if(auto *pContext = m_pPlot->context(); !pContext) return;
		QPainter painter(m_pPlot);
		painter.setRenderHint(QPainter::Antialiasing);

		// gl painting
		{
			BOOST_SCOPE_EXIT(&painter) { painter.endNativePainting(); } BOOST_SCOPE_EXIT_END

			if(m_bPickerNeedsUpdate) UpdatePicker();

			auto *pGl = GetGlFunctions();
			painter.beginNativePainting();
			DoPaintGL(pGl);
		}

		// qt painting
		DoPaintNonGL(painter);
	}
	else	// threaded
	{
		QThread *pThisThread = QThread::currentThread();
		if(!pThisThread->isRunning() || pThisThread->isInterruptionRequested())
			return;

		if(auto *pContext = m_pPlot->context(); !pContext) return;

#if QT_VERSION >= QT_VERSION_CHECK(5, 10, 0)
		QMetaObject::invokeMethod(m_pPlot, &GlPlot::MoveContextToThread, Qt::ConnectionType::BlockingQueuedConnection);
#else
		QMetaObject::invokeMethod(m_pPlot, "MoveContextToThread", Qt::ConnectionType::BlockingQueuedConnection);
#endif
		if(!m_pPlot->IsContextInThread())
		{
			std::cerr << __func__ << ": Context is not in thread!" << std::endl;
			return;
		}

		m_pPlot->GetMutex()->lock();
		m_pPlot->makeCurrent();
		BOOST_SCOPE_EXIT(m_pPlot)
		{
			m_pPlot->doneCurrent();
			m_pPlot->context()->moveToThread(qGuiApp->thread());
			m_pPlot->GetMutex()->unlock();

			if constexpr(!m_usetimer)
			{
				// if the frame is not already updated by the timer, directly update it
				m_pPlot->GetRenderer()->RequestPlotUpdate();
			}
		}
		BOOST_SCOPE_EXIT_END

		if(!m_bInitialised) initialiseGL();
		if(!m_bInitialised)
		{
			std::cerr << "Cannot initialise GL." << std::endl;
			return;
		}

		if(m_bWantsResize) resizeGL();
		if(m_bPickerNeedsUpdate) UpdatePicker();

		DoPaintGL(GetGlFunctions());
	}
}
// ----------------------------------------------------------------------------




// ----------------------------------------------------------------------------
// GLPlot wrapper class
// ----------------------------------------------------------------------------

GlPlot::GlPlot(QWidget *pParent) : QOpenGLWidget(pParent),
	m_renderer(std::make_unique<GlPlotRenderer>(this)),
	m_thread_impl(std::make_unique<QThread>(this))
{
	qRegisterMetaType<std::size_t>("std::size_t");

	if constexpr(m_isthreaded)
	{
		m_renderer->moveToThread(m_thread_impl.get());

		connect(m_thread_impl.get(), &QThread::started, m_renderer.get(), &GlPlotRenderer::startedThread);
		connect(m_thread_impl.get(), &QThread::finished, m_renderer.get(), &GlPlotRenderer::stoppedThread);
	}

	connect(this, &QOpenGLWidget::aboutToCompose, this, &GlPlot::beforeComposing);
	connect(this, &QOpenGLWidget::frameSwapped, this, &GlPlot::afterComposing);
	connect(this, &QOpenGLWidget::aboutToResize, this, &GlPlot::beforeResizing);
	connect(this, &QOpenGLWidget::resized, this, &GlPlot::afterResizing);

	//setUpdateBehavior(QOpenGLWidget::PartialUpdate);
	setMouseTracking(true);

	if constexpr(m_isthreaded)
		m_thread_impl->start();
}


GlPlot::~GlPlot()
{
	setMouseTracking(false);

	if constexpr(m_isthreaded)
	{
		m_thread_impl->requestInterruption();
		m_thread_impl->exit();
		m_thread_impl->wait();
	}
}


void GlPlot::initializeGL()
{
	if constexpr(!m_isthreaded)
	{
		m_renderer->initialiseGL();
		if(m_renderer->IsInitialised())
			emit AfterGLInitialisation();
		else
			emit GLInitialisationFailed();
	}
}


void GlPlot::resizeGL(int w, int h)
{
	if constexpr(!m_isthreaded)
	{
		m_renderer->SetScreenDims(w, h);
		m_renderer->resizeGL();
	}
}


void GlPlot::paintGL()
{
	if constexpr(!m_isthreaded)
	{
		m_renderer->paintGL();
	}
}


void GlPlot::mouseMoveEvent(QMouseEvent *pEvt)
{
	m_renderer->mouseMoveEvent(pEvt->localPos());
	m_mouseMovedBetweenDownAndUp = 1;
	pEvt->accept();
}


void GlPlot::mousePressEvent(QMouseEvent *pEvt)
{
	m_mouseMovedBetweenDownAndUp = 0;

	if(pEvt->buttons() & Qt::LeftButton) m_mouseDown[0] = 1;
	if(pEvt->buttons() & Qt::MiddleButton) m_mouseDown[1] = 1;
	if(pEvt->buttons() & Qt::RightButton) m_mouseDown[2] = 1;

	if(m_mouseDown[1])
		m_renderer->ResetZoom();
	if(m_mouseDown[2])
		m_renderer->BeginRotation();

	pEvt->accept();
	emit MouseDown(m_mouseDown[0], m_mouseDown[1], m_mouseDown[2]);
}


void GlPlot::mouseReleaseEvent(QMouseEvent *pEvt)
{
	bool mouseDownOld[] = { m_mouseDown[0], m_mouseDown[1], m_mouseDown[2] };

	if((pEvt->buttons() & Qt::LeftButton) == 0) m_mouseDown[0] = 0;
	if((pEvt->buttons() & Qt::MiddleButton) == 0) m_mouseDown[1] = 0;
	if((pEvt->buttons() & Qt::RightButton) == 0) m_mouseDown[2] = 0;

	if(!m_mouseDown[2])
		m_renderer->EndRotation();

	pEvt->accept();
	emit MouseUp(!m_mouseDown[0], !m_mouseDown[1], !m_mouseDown[2]);

	// only emit click if moving the mouse (i.e. rotationg the scene) was not the primary intent
	if(!m_mouseMovedBetweenDownAndUp)
	{
		bool mouseClicked[] = { !m_mouseDown[0] && mouseDownOld[0],
			!m_mouseDown[1] && mouseDownOld[1],
			!m_mouseDown[2] && mouseDownOld[2] };
		if(mouseClicked[0] || mouseClicked[1] || mouseClicked[2])
			emit MouseClick(mouseClicked[0], mouseClicked[1], mouseClicked[2]);
	}
}


void GlPlot::wheelEvent(QWheelEvent *pEvt)
{
	const t_real_gl degrees = pEvt->angleDelta().y() / 8.;
	m_renderer->zoom(degrees);

	pEvt->accept();
}


void GlPlot::paintEvent(QPaintEvent* pEvt)
{
	if constexpr(!m_isthreaded)
		QOpenGLWidget::paintEvent(pEvt);
}

/**
 * move the GL context to the associated thread
 */
void GlPlot::MoveContextToThread()
{
	if constexpr(m_isthreaded)
	{
		if(auto *pContext = context(); pContext && m_thread_impl.get())
			pContext->moveToThread(m_thread_impl.get());
	}
}


/**
 * does the GL context run in the current thread?
 */
bool GlPlot::IsContextInThread() const
{
	if constexpr(m_isthreaded)
	{
		auto *pContext = context();
		if(!pContext) return false;

		return pContext->thread() == m_thread_impl.get();
	}
	else
	{
		return true;
	}
}


/**
 * main thread wants to compose -> wait for sub-threads to be finished
 */
void GlPlot::beforeComposing()
{
	if constexpr(m_isthreaded)
	{
		m_mutex.lock();
	}
}


/**
 * main thread has composed -> sub-threads can be unblocked
 */
void GlPlot::afterComposing()
{
	if constexpr(m_isthreaded)
	{
		m_mutex.unlock();
#if QT_VERSION >= QT_VERSION_CHECK(5, 10, 0)
		QMetaObject::invokeMethod(m_renderer.get(), &GlPlotRenderer::paintGL, Qt::ConnectionType::QueuedConnection);
#else
		QMetaObject::invokeMethod(m_renderer.get(), "paintGL", Qt::ConnectionType::QueuedConnection);
#endif
	}
}


/**
 * main thread wants to resize -> wait for sub-threads to be finished
 */
void GlPlot::beforeResizing()
{
	if constexpr(m_isthreaded)
	{
		m_mutex.lock();
	}
}


/**
 * main thread has resized -> sub-threads can be unblocked
 */
void GlPlot::afterResizing()
{
	if constexpr(m_isthreaded)
	{
		m_mutex.unlock();

		const int w = width(), h = height();
		m_renderer->SetScreenDims(w, h);
	}
}

// ----------------------------------------------------------------------------
}
