/**
 * GL plotter
 * @author Tobias Weber <tweber@ill.fr>
 * @date Nov-2017 -- 2019
 * @license GPLv3, see 'LICENSE' file
 * @desc The present version was forked on 8-Nov-2018 from the privately developed "magtools" project (https://github.com/t-weber/magtools).
 */

#include "glplot.h"

#include <QtGui/QOpenGLContext>
#include <QtGui/QOpenGLFunctions>
#include <QtGui/QSurfaceFormat>
#include <QtGui/QPainter>
#include <QtGui/QGuiApplication>

#include <iostream>
#include <boost/scope_exit.hpp>
#include <boost/algorithm/string/replace.hpp>
namespace algo = boost::algorithm;



// ----------------------------------------------------------------------------
void set_gl_format(bool bCore, int iMajorVer, int iMinorVer, int iSamples)
{
	QSurfaceFormat surf = QSurfaceFormat::defaultFormat();

	surf.setRenderableType(QSurfaceFormat::OpenGL);
	if(bCore)
		surf.setProfile(QSurfaceFormat::CoreProfile);
	else
		surf.setProfile(QSurfaceFormat::CompatibilityProfile);

	if(iMajorVer > 0 && iMinorVer > 0)
		surf.setVersion(iMajorVer, iMinorVer);

	surf.setSwapBehavior(QSurfaceFormat::DoubleBuffer);
	if(iSamples > 0)
		surf.setSamples(iSamples);	// multisampling

	QSurfaceFormat::setDefaultFormat(surf);
}
// ----------------------------------------------------------------------------




// ----------------------------------------------------------------------------
// GL plot implementation

GlPlot_impl::GlPlot_impl(GlPlot *pPlot) : m_pPlot{pPlot}
{
	if constexpr(m_usetimer)
	{
		connect(&m_timer, &QTimer::timeout,
			this, static_cast<void (GlPlot_impl::*)()>(&GlPlot_impl::tick));
		m_timer.start(std::chrono::milliseconds(1000 / 60));
	}

	UpdateCam();
}


GlPlot_impl::~GlPlot_impl()
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
	{
		obj.m_pvertexbuf.reset();
		obj.m_pnormalsbuf.reset();
		obj.m_pcolorbuf.reset();
		if(pGl) pGl->glDeleteVertexArrays(1, &obj.m_vertexarr);
	}

	m_objs.clear();
	LOGGLERR(pGl)
}


void GlPlot_impl::startedThread() { }
void GlPlot_impl::stoppedThread() { }


qgl_funcs* GlPlot_impl::GetGlFunctions(QOpenGLWidget *pWidget)
{
	if(!pWidget) pWidget = (QOpenGLWidget*)m_pPlot;
	qgl_funcs *pGl = nullptr;

	if constexpr(std::is_same_v<qgl_funcs, QOpenGLFunctions>)
		pGl = (qgl_funcs*)pWidget->context()->functions();
	else
		pGl = (qgl_funcs*)pWidget->context()->versionFunctions<qgl_funcs>();

	if(!pGl)
		std::cerr << "No suitable GL interface found." << std::endl;

	return pGl;
}


QPointF GlPlot_impl::GlToScreenCoords(const t_vec_gl& vec4, bool *pVisible)
{
	auto [ vecPersp, vec ] =
	m::hom_to_screen_coords<t_mat_gl, t_vec_gl>(vec4, m_matCam, m_matPerspective, m_matViewport, true);

	// position not visible -> return a point outside the viewport
	if(vecPersp[2] > 1.)
	{
		if(pVisible) *pVisible = false;
		return QPointF(-1*m_iScreenDims[0], -1*m_iScreenDims[1]);
	}

	if(pVisible) *pVisible = true;
	return QPointF(vec[0], vec[1]);
}


t_mat_gl GlPlot_impl::GetArrowMatrix(const t_vec_gl& vecTo, t_real_gl postscale, const t_vec_gl& vecPostTrans,
	const t_vec_gl& vecFrom, t_real_gl prescale, const t_vec_gl& vecPreTrans)
{
	t_mat_gl mat = m::unit<t_mat_gl>(4);

	mat *= m::hom_translation<t_mat_gl>(vecPreTrans[0], vecPreTrans[1], vecPreTrans[2]);
	mat *= m::hom_scaling<t_mat_gl>(prescale, prescale, prescale);

	mat *= m::rotation<t_mat_gl, t_vec_gl>(vecFrom, vecTo);

	mat *= m::hom_scaling<t_mat_gl>(postscale, postscale, postscale);
	mat *= m::hom_translation<t_mat_gl>(vecPostTrans[0], vecPostTrans[1], vecPostTrans[2]);

	return mat;
}


GlPlotObj GlPlot_impl::CreateTriangleObject(const std::vector<t_vec3_gl>& verts,
	const std::vector<t_vec3_gl>& triagverts, const std::vector<t_vec3_gl>& norms,
	const t_vec_gl& color, bool bUseVertsAsNorm)
{
	// TODO: move context to calling thread
	m_pPlot->makeCurrent();
	BOOST_SCOPE_EXIT(m_pPlot) { m_pPlot->doneCurrent(); } BOOST_SCOPE_EXIT_END


	qgl_funcs* pGl = GetGlFunctions();

	GLint attrVertex = m_attrVertex;
	GLint attrVertexNormal = m_attrVertexNorm;
	GLint attrVertexColor = m_attrVertexCol;

	GlPlotObj obj;
	obj.m_type = GlPlotObjType::TRIANGLES;
	obj.m_color = color;

	// flatten vertex array into raw float array
	auto to_float_array = [](const std::vector<t_vec3_gl>& verts, int iRepeat=1, int iElems=3, bool bNorm=false)
		-> std::vector<t_real_gl>
	{
		std::vector<t_real_gl> vecRet;
		vecRet.reserve(iRepeat*verts.size()*iElems);

		for(const t_vec3_gl& vert : verts)
		{
			t_real_gl norm = bNorm ? m::norm<t_vec3_gl>(vert) : 1;

			for(int i=0; i<iRepeat; ++i)
				for(int iElem=0; iElem<iElems; ++iElem)
					vecRet.push_back(vert[iElem] / norm);
		}

		return vecRet;
	};

	// main vertex array object
	pGl->glGenVertexArrays(1, &obj.m_vertexarr);
	pGl->glBindVertexArray(obj.m_vertexarr);

	{	// vertices
		obj.m_pvertexbuf = std::make_shared<QOpenGLBuffer>(QOpenGLBuffer::VertexBuffer);

		if(!obj.m_pvertexbuf->create())
			std::cerr << "Cannot create vertex buffer." << std::endl;
		if(!obj.m_pvertexbuf->bind())
			std::cerr << "Cannot bind vertex buffer." << std::endl;
		BOOST_SCOPE_EXIT(&obj) { obj.m_pvertexbuf->release(); } BOOST_SCOPE_EXIT_END

		auto vecVerts = to_float_array(triagverts, 1,3, false);
		obj.m_pvertexbuf->allocate(vecVerts.data(), vecVerts.size()*sizeof(typename decltype(vecVerts)::value_type));
		pGl->glVertexAttribPointer(attrVertex, 3, GL_FLOAT, 0, 0, nullptr);
	}

	{	// normals
		obj.m_pnormalsbuf = std::make_shared<QOpenGLBuffer>(QOpenGLBuffer::VertexBuffer);

		obj.m_pnormalsbuf->create();
		obj.m_pnormalsbuf->bind();
		BOOST_SCOPE_EXIT(&obj) { obj.m_pnormalsbuf->release(); } BOOST_SCOPE_EXIT_END

		auto vecNorms = bUseVertsAsNorm ? to_float_array(triagverts, 1,3, true) : to_float_array(norms, 3,3, false);
		obj.m_pnormalsbuf->allocate(vecNorms.data(), vecNorms.size()*sizeof(typename decltype(vecNorms)::value_type));
		pGl->glVertexAttribPointer(attrVertexNormal, 3, GL_FLOAT, 0, 0, nullptr);
	}

	{	// colors
		obj.m_pcolorbuf = std::make_shared<QOpenGLBuffer>(QOpenGLBuffer::VertexBuffer);

		obj.m_pcolorbuf->create();
		obj.m_pcolorbuf->bind();
		BOOST_SCOPE_EXIT(&obj) { obj.m_pcolorbuf->release(); } BOOST_SCOPE_EXIT_END

		std::vector<t_real_gl> vecCols;
		vecCols.reserve(4*triagverts.size());
		for(std::size_t iVert=0; iVert<triagverts.size(); ++iVert)
		{
			for(int icol=0; icol<obj.m_color.size(); ++icol)
				vecCols.push_back(obj.m_color[icol]);
		}

		obj.m_pcolorbuf->allocate(vecCols.data(), vecCols.size()*sizeof(typename decltype(vecCols)::value_type));
		pGl->glVertexAttribPointer(attrVertexColor, 4, GL_FLOAT, 0, 0, nullptr);
	}


	obj.m_vertices = std::move(verts);
	obj.m_triangles = std::move(triagverts);
	LOGGLERR(pGl)

	return obj;
}


GlPlotObj GlPlot_impl::CreateLineObject(const std::vector<t_vec3_gl>& verts, const t_vec_gl& color)
{
	// TODO: move context to calling thread
	m_pPlot->makeCurrent();
	BOOST_SCOPE_EXIT(m_pPlot) { m_pPlot->doneCurrent(); } BOOST_SCOPE_EXIT_END


	qgl_funcs* pGl = GetGlFunctions();
	GLint attrVertex = m_attrVertex;
	GLint attrVertexColor = m_attrVertexCol;

	GlPlotObj obj;
	obj.m_type = GlPlotObjType::LINES;
	obj.m_color = color;

	// flatten vertex array into raw float array
	auto to_float_array = [](const std::vector<t_vec3_gl>& verts, int iElems=3) -> std::vector<t_real_gl>
	{
		std::vector<t_real_gl> vecRet;
		vecRet.reserve(verts.size()*iElems);

		for(const t_vec3_gl& vert : verts)
		{
			for(int iElem=0; iElem<iElems; ++iElem)
				vecRet.push_back(vert[iElem]);
		}

		return vecRet;
	};

	// main vertex array object
	pGl->glGenVertexArrays(1, &obj.m_vertexarr);
	pGl->glBindVertexArray(obj.m_vertexarr);

	{	// vertices
		obj.m_pvertexbuf = std::make_shared<QOpenGLBuffer>(QOpenGLBuffer::VertexBuffer);

		obj.m_pvertexbuf->create();
		obj.m_pvertexbuf->bind();
		BOOST_SCOPE_EXIT(&obj) { obj.m_pvertexbuf->release(); } BOOST_SCOPE_EXIT_END

		auto vecVerts = to_float_array(verts, 3);
		obj.m_pvertexbuf->allocate(vecVerts.data(), vecVerts.size()*sizeof(typename decltype(vecVerts)::value_type));
		pGl->glVertexAttribPointer(attrVertex, 3, GL_FLOAT, 0, 0, nullptr);
	}

	{	// colors
		obj.m_pcolorbuf = std::make_shared<QOpenGLBuffer>(QOpenGLBuffer::VertexBuffer);

		obj.m_pcolorbuf->create();
		obj.m_pcolorbuf->bind();
		BOOST_SCOPE_EXIT(&obj) { obj.m_pcolorbuf->release(); } BOOST_SCOPE_EXIT_END

		std::vector<t_real_gl> vecCols;
		vecCols.reserve(4*verts.size());
		for(std::size_t iVert=0; iVert<verts.size(); ++iVert)
		{
			for(int icol=0; icol<obj.m_color.size(); ++icol)
				vecCols.push_back(obj.m_color[icol]);
		}

		obj.m_pcolorbuf->allocate(vecCols.data(), vecCols.size()*sizeof(typename decltype(vecCols)::value_type));
		pGl->glVertexAttribPointer(attrVertexColor, 4, GL_FLOAT, 0, 0, nullptr);
	}


	obj.m_vertices = std::move(verts);
	LOGGLERR(pGl)

	return obj;
}


void GlPlot_impl::SetObjectMatrix(std::size_t idx, const t_mat_gl& mat)
{
	if(idx >= m_objs.size()) return;
	m_objs[idx].m_mat = mat;
}


void GlPlot_impl::SetObjectCol(std::size_t idx, t_real_gl r, t_real_gl g, t_real_gl b, t_real_gl a)
{
	if(idx >= m_objs.size()) return;
	m_objs[idx].m_color = m::create<t_vec_gl>({r,g,b,a});
}

void GlPlot_impl::SetObjectLabel(std::size_t idx, const std::string& label)
{
	if(idx >= m_objs.size()) return;
	m_objs[idx].m_label = label;
}


void GlPlot_impl::SetObjectDataString(std::size_t idx, const std::string& data)
{
	if(idx >= m_objs.size()) return;
	m_objs[idx].m_datastr = data;
}

const std::string& GlPlot_impl::GetObjectDataString(std::size_t idx) const
{
	static const std::string empty{};
	if(idx >= m_objs.size()) return empty;

	return m_objs[idx].m_datastr;
}


void GlPlot_impl::SetObjectVisible(std::size_t idx, bool visible)
{
	if(idx >= m_objs.size()) return;
	m_objs[idx].m_visible = visible;
}


void GlPlot_impl::SetObjectHighlight(std::size_t idx, bool highlight)
{
	if(idx >= m_objs.size()) return;
	m_objs[idx].m_highlighted = highlight;
}


bool GlPlot_impl::GetObjectHighlight(std::size_t idx) const
{
	if(idx >= m_objs.size()) return 0;
	return m_objs[idx].m_highlighted;
}


void GlPlot_impl::RemoveObject(std::size_t idx)
{
	m_objs[idx].m_valid = false;

	// TODO: remove if object has no follow-up indices
}


std::size_t GlPlot_impl::AddLinkedObject(std::size_t linkTo,
	t_real_gl x, t_real_gl y, t_real_gl z,
	t_real_gl r, t_real_gl g, t_real_gl b, t_real_gl a)
{
	GlPlotObj obj;
	obj.linkedObj = linkTo;
	obj.m_mat = m::hom_translation<t_mat_gl>(x, y, z);
	obj.m_color = m::create<t_vec_gl>({r, g, b, a});

	QMutexLocker _locker{&m_mutexObj};
	m_objs.emplace_back(std::move(obj));

	return m_objs.size()-1;		// object handle
}


std::size_t GlPlot_impl::AddSphere(t_real_gl rad, t_real_gl x, t_real_gl y, t_real_gl z,
	t_real_gl r, t_real_gl g, t_real_gl b, t_real_gl a)
{
	constexpr int numsubdivs = 2;

	auto solid = m::create_icosahedron<t_vec3_gl>(1);
	auto [triagverts, norms, uvs] = m::spherify<t_vec3_gl>(
		m::subdivide_triangles<t_vec3_gl>(
			m::create_triangles<t_vec3_gl>(solid), numsubdivs), rad);
	auto [boundingSpherePos, boundingSphereRad] = m::bounding_sphere<t_vec3_gl>(triagverts);

	QMutexLocker _locker{&m_mutexObj};

	auto obj = CreateTriangleObject(std::get<0>(solid), triagverts, norms, m::create<t_vec_gl>({r,g,b,a}), true);
	obj.m_mat = m::hom_translation<t_mat_gl>(x, y, z);
	obj.m_boundingSpherePos = std::move(boundingSpherePos);
	obj.m_boundingSphereRad = boundingSphereRad;
	//obj.m_boundingSphereRad = rad;
	m_objs.emplace_back(std::move(obj));

	return m_objs.size()-1;		// object handle
}


std::size_t GlPlot_impl::AddCylinder(t_real_gl rad, t_real_gl h,
	t_real_gl x, t_real_gl y, t_real_gl z,
	t_real_gl r, t_real_gl g, t_real_gl b, t_real_gl a)
{
	auto solid = m::create_cylinder<t_vec3_gl>(rad, h, true);
	auto [triagverts, norms, uvs] = m::create_triangles<t_vec3_gl>(solid);
	auto [boundingSpherePos, boundingSphereRad] = m::bounding_sphere<t_vec3_gl>(triagverts);

	QMutexLocker _locker{&m_mutexObj};

	auto obj = CreateTriangleObject(std::get<0>(solid), triagverts, norms, m::create<t_vec_gl>({r,g,b,a}), false);
	obj.m_mat = m::hom_translation<t_mat_gl>(x, y, z);
	obj.m_boundingSpherePos = std::move(boundingSpherePos);
	obj.m_boundingSphereRad = boundingSphereRad;
	m_objs.emplace_back(std::move(obj));

	return m_objs.size()-1;		// object handle
}


std::size_t GlPlot_impl::AddCone(t_real_gl rad, t_real_gl h,
	t_real_gl x, t_real_gl y, t_real_gl z,
	t_real_gl r, t_real_gl g, t_real_gl b, t_real_gl a)
{
	auto solid = m::create_cone<t_vec3_gl>(rad, h);
	auto [triagverts, norms, uvs] = m::create_triangles<t_vec3_gl>(solid);
	auto [boundingSpherePos, boundingSphereRad] = m::bounding_sphere<t_vec3_gl>(triagverts);

	QMutexLocker _locker{&m_mutexObj};

	auto obj = CreateTriangleObject(std::get<0>(solid), triagverts, norms, m::create<t_vec_gl>({r,g,b,a}), false);
	obj.m_mat = m::hom_translation<t_mat_gl>(x, y, z);
	obj.m_boundingSpherePos = std::move(boundingSpherePos);
	obj.m_boundingSphereRad = boundingSphereRad;
	m_objs.emplace_back(std::move(obj));

	return m_objs.size()-1;		// object handle
}


std::size_t GlPlot_impl::AddArrow(t_real_gl rad, t_real_gl h,
	t_real_gl x, t_real_gl y, t_real_gl z,
	t_real_gl r, t_real_gl g, t_real_gl b, t_real_gl a)
{
	auto solid = m::create_cylinder<t_vec3_gl>(rad, h, 2, 32, rad, rad*1.5);
	auto [triagverts, norms, uvs] = m::create_triangles<t_vec3_gl>(solid);
	auto [boundingSpherePos, boundingSphereRad] = m::bounding_sphere<t_vec3_gl>(triagverts);

	QMutexLocker _locker{&m_mutexObj};

	auto obj = CreateTriangleObject(std::get<0>(solid), triagverts, norms, m::create<t_vec_gl>({r,g,b,a}), false);
	obj.m_mat = GetArrowMatrix(m::create<t_vec_gl>({1,0,0}), 1., m::create<t_vec_gl>({x,y,z}), m::create<t_vec_gl>({0,0,1}));
	obj.m_boundingSpherePos = std::move(boundingSpherePos);
	obj.m_boundingSphereRad = boundingSphereRad;
	obj.m_labelPos = m::create<t_vec3_gl>({0., 0., 0.75});
	m_objs.emplace_back(std::move(obj));

	return m_objs.size()-1;		// object handle
}


std::size_t GlPlot_impl::AddCoordinateCross(t_real_gl min, t_real_gl max)
{
	auto col = m::create<t_vec_gl>({0,0,0,1});
	auto verts = std::vector<t_vec3_gl>
	{{
		m::create<t_vec3_gl>({min,0,0}), m::create<t_vec3_gl>({max,0,0}),
		m::create<t_vec3_gl>({0,min,0}), m::create<t_vec3_gl>({0,max,0}),
		m::create<t_vec3_gl>({0,0,min}), m::create<t_vec3_gl>({0,0,max}),
	}};

	QMutexLocker _locker{&m_mutexObj};

	auto obj = CreateLineObject(verts, col);
	obj.m_invariant = true;
	m_objs.emplace_back(std::move(obj));

	return m_objs.size()-1;		// object handle
}



void GlPlot_impl::initialiseGL()
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
 * phong lighting model, see: https://en.wikipedia.org/wiki/Phong_reflection_model
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
	fragcol = vertexcol * I;
	fragcol *= constcol;
	fragcol[3] = 1;
})RAW";
// --------------------------------------------------------------------


	// set glsl version and constants
	const std::string strGlsl = std::to_string(_GL_MAJ_VER*100 + _GL_MIN_VER*10);
	const std::string strPi = std::to_string(m::pi<t_real_gl>);
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


	// options
	pGl->glCullFace(GL_BACK);
	pGl->glEnable(GL_CULL_FACE);

	pGl->glEnable(GL_MULTISAMPLE);
	pGl->glEnable(GL_LINE_SMOOTH);
	pGl->glEnable(GL_POLYGON_SMOOTH);
	pGl->glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
	pGl->glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);

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


void GlPlot_impl::SetScreenDims(int w, int h)
{
	m_iScreenDims[0] = w;
	m_iScreenDims[1] = h;
	m_bWantsResize = true;
}


void GlPlot_impl::resizeGL()
{
	if(!m_bPlatformSupported || !m_bInitialised) return;

	const int w = m_iScreenDims[0];
	const int h = m_iScreenDims[1];

	if(auto *pContext = ((QOpenGLWidget*)m_pPlot)->context(); !pContext)
		return;
	auto *pGl = GetGlFunctions();
	if(!pGl)
		return;

	m_matViewport = m::hom_viewport<t_mat_gl>(w, h, 0., 1.);
	std::tie(m_matViewport_inv, std::ignore) = m::inv<t_mat_gl>(m_matViewport);

	m_matPerspective = m::hom_perspective<t_mat_gl>(0.01, 100., m::pi<t_real_gl>*0.5, t_real_gl(h)/t_real_gl(w));
	//m_matPerspective = m::hom_ortho<t_mat_gl>(0.01, 100., -t_real_gl(w)*0.0025, t_real_gl(w)*0.0025, -t_real_gl(h)*0.0025, t_real_gl(h)*0.0025);
	std::tie(m_matPerspective_inv, std::ignore) = m::inv<t_mat_gl>(m_matPerspective);

	pGl->glViewport(0, 0, w, h);
	pGl->glDepthRange(0, 1);

	// bind shaders
	m_pShaders->bind();
	BOOST_SCOPE_EXIT(m_pShaders) { m_pShaders->release(); } BOOST_SCOPE_EXIT_END
	LOGGLERR(pGl);

	// set matrices
	m_pShaders->setUniformValue(m_uniMatrixCam, m_matCam);
	m_pShaders->setUniformValue(m_uniMatrixCamInv, m_matCam_inv);
	m_pShaders->setUniformValue(m_uniMatrixProj, m_matPerspective);
	LOGGLERR(pGl);

	m_bWantsResize = false;
}


/**
 * set up a (crystal) B matrix
 */
void GlPlot_impl::SetBTrafo(const t_mat_gl& matB, const t_mat_gl* matA)
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
		std::tie(m_matA, ok) = m::inv(m_matB);
		if(!ok)
		{
			m_matA = m::unit<t_mat_gl>();
			std::cerr << "Error: Cannot invert B matrix." << std::endl;
		}
		else
		{
			m_matA *= t_real_gl(2)*m::pi<t_real_gl>;
			m_matA(3,3) = 1;
		}
	}

	m_bBTrafoNeedsUpdate = true;
	RequestPlotUpdate();
}


void GlPlot_impl::SetCoordSys(int iSys)
{
	m_iCoordSys = iSys;
	RequestPlotUpdate();
}


/**
 * update the shader's B matrix
 */
void GlPlot_impl::UpdateBTrafo()
{
	m_pShaders->setUniformValue(m_uniMatrixA, m_matA);
	m_pShaders->setUniformValue(m_uniMatrixB, m_matB);

	m_bBTrafoNeedsUpdate = false;
}


void GlPlot_impl::UpdateCam()
{
	m_matCam = m_matCamBase;
	m_matCam(2,3) /= m_zoom;
	m_matCam *= m_matCamRot;
	std::tie(m_matCam_inv, std::ignore) = m::inv<t_mat_gl>(m_matCam);

	/*auto M = m_matCam_inv; M(0,3) = M(1,3) = M(2,3) = 0;
	t_vec3_gl vecCamPos = M * m::create<t_vec3_gl>({-m_matCam(0,3), -m_matCam(1,3), -m_matCam(2,3)});
	std::cout << vecCamPos[0] << " " << vecCamPos[1] << " " << vecCamPos[2] << std::endl;*/

	m_bPickerNeedsUpdate = true;
	RequestPlotUpdate();
}


/**
 * request a plot update
 */
void GlPlot_impl::RequestPlotUpdate()
{
	QMetaObject::invokeMethod((QOpenGLWidget*)m_pPlot,
		static_cast<void (QOpenGLWidget::*)()>(&QOpenGLWidget::update),
		Qt::ConnectionType::QueuedConnection);
}


void GlPlot_impl::SetLight(std::size_t idx, const t_vec3_gl& pos)
{
	if(m_lights.size() < idx+1)
		m_lights.resize(idx+1);

	m_lights[idx] = pos;
	m_bLightsNeedUpdate = true;
}


void GlPlot_impl::UpdateLights()
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


void GlPlot_impl::EnablePicker(bool b)
{
	m_bPickerEnabled = b;
}


void GlPlot_impl::UpdatePicker()
{
	if(!m_bInitialised || !m_bPlatformSupported || !m_bPickerEnabled) return;

	// picker ray
	auto [org, dir] = m::hom_line_from_screen_coords<t_mat_gl, t_vec_gl>(
		m_posMouse.x(), m_posMouse.y(), 0., 1., m_matCam_inv,
		m_matPerspective_inv, m_matViewport_inv, &m_matViewport, true);
	t_vec3_gl org3 = m::create<t_vec3_gl>({org[0], org[1], org[2]});
	t_vec3_gl dir3 = m::create<t_vec3_gl>({dir[0], dir[1], dir[2]});


	// intersection with unit sphere around origin
	bool hasSphereInters = false;
	t_vec_gl vecClosestSphereInters = m::create<t_vec_gl>({0,0,0,0});

	auto intersUnitSphere =
	m::intersect_line_sphere<t_vec3_gl, std::vector>(org3, dir3,
		m::create<t_vec3_gl>({0,0,0}), t_real_gl(m_pickerSphereRadius));
	for(const auto& result : intersUnitSphere)
	{
		t_vec_gl vecInters4 = m::create<t_vec_gl>({result[0], result[1], result[2], 1});

		if(!hasSphereInters)
		{	// first intersection
			vecClosestSphereInters = vecInters4;
			hasSphereInters = true;
		}
		else
		{	// test if next intersection is closer...
			t_vec_gl oldPosTrafo = m_matCam * vecClosestSphereInters;
			t_vec_gl newPosTrafo = m_matCam * vecInters4;

			// ... it is closer.
			if(m::norm(newPosTrafo) < m::norm(oldPosTrafo))
				vecClosestSphereInters = vecInters4;
		}
	}


	// crystal or lab coordinate system?
	const t_mat_gl matUnit = m::unit<t_mat_gl>();
	const t_mat_gl *coordTrafo = &matUnit;
	t_mat_gl coordTrafoInv = matUnit;
	if(m_iCoordSys == 1)
	{
		coordTrafo = &m_matA;
		coordTrafoInv = m_matB / (t_real_gl(2)*m::pi<t_real_gl>);
		coordTrafoInv(3,3) = 1;
	}


	// intersection with geometry
	bool hasInters = false;
	t_vec_gl vecClosestInters = m::create<t_vec_gl>({0,0,0,0});
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


		t_mat_gl matTrafo = (*coordTrafo) * obj.m_mat * coordTrafoInv;

		// scaling factor, TODO: maximum factor for non-uniform scaling
		auto scale = std::cbrt(std::abs(m::det(matTrafo)));

		// intersection with bounding sphere?
		auto boundingInters =
			m::intersect_line_sphere<t_vec3_gl, std::vector>(org3, dir3,
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

			// coordTrafoInv only keeps 3d objects from locally distorting
			auto [vecInters, bInters, lamInters] =
				m::intersect_line_poly<t_vec3_gl, t_mat_gl>(org3, dir3, poly, matTrafo);

			if(bInters)
			{
				t_vec_gl vecInters4 = m::create<t_vec_gl>({vecInters[0], vecInters[1], vecInters[2], 1});

				if(!hasInters)
				{	// first intersection
					vecClosestInters = vecInters4;
					objInters = curObj;
					hasInters = true;
				}
				else
				{	// test if next intersection is closer...
					t_vec_gl oldPosTrafo = m_matCam * vecClosestInters;
					t_vec_gl newPosTrafo = m_matCam * vecInters4;

					if(m::norm(newPosTrafo) < m::norm(oldPosTrafo))
					{	// ...it is closer
						vecClosestInters = vecInters4;
						objInters = curObj;
					}
				}
			}
		}
	}

	m_bPickerNeedsUpdate = false;
	t_vec3_gl vecClosestInters3 = m::create<t_vec3_gl>({vecClosestInters[0], vecClosestInters[1], vecClosestInters[2]});
	t_vec3_gl vecClosestSphereInters3 = m::create<t_vec3_gl>({vecClosestSphereInters[0], vecClosestSphereInters[1], vecClosestSphereInters[2]});
	emit PickerIntersection(hasInters ? &vecClosestInters3 : nullptr, objInters,
		hasSphereInters ? &vecClosestSphereInters3 : nullptr);
}


void GlPlot_impl::mouseMoveEvent(const QPointF& pos)
{
	m_posMouse = pos;

	if(m_bInRotation)
	{
		auto diff = (m_posMouse - m_posMouseRotationStart);
		t_real_gl phi = diff.x() + m_phi_saved;
		t_real_gl theta = diff.y() + m_theta_saved;

		m_matCamRot = m::rotation<t_mat_gl, t_vec_gl>(m_vecCamX, theta/180.*M_PI, 0);
		m_matCamRot *= m::rotation<t_mat_gl, t_vec_gl>(m_vecCamY, phi/180.*M_PI, 0);

		UpdateCam();
	}
	else
	{
		// also automatically done in UpdateCam
		m_bPickerNeedsUpdate = true;
		RequestPlotUpdate();
	}
}


void GlPlot_impl::zoom(t_real_gl val)
{
	m_zoom *= std::pow(2., val/64.);
	UpdateCam();
}


void GlPlot_impl::ResetZoom()
{
	m_zoom = 1;
	UpdateCam();
}


void GlPlot_impl::BeginRotation()
{
	if(!m_bInRotation)
	{
		m_posMouseRotationStart = m_posMouse;
		m_bInRotation = true;
	}
}


void GlPlot_impl::EndRotation()
{
	if(m_bInRotation)
	{
		auto diff = (m_posMouse - m_posMouseRotationStart);
		m_phi_saved += diff.x();
		m_theta_saved += diff.y();

		m_bInRotation = false;
	}
}


void GlPlot_impl::tick()
{
	tick(std::chrono::milliseconds(1000 / 60));
}


void GlPlot_impl::tick(const std::chrono::milliseconds& ms)
{
	// TODO
	UpdateCam();
}


/**
 * pure gl drawing
 */
void GlPlot_impl::DoPaintGL(qgl_funcs *pGl)
{
	if(!pGl)
		return;

	// clear
	pGl->glClearColor(1., 1., 1., 1.);
	pGl->glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	pGl->glEnable(GL_DEPTH_TEST);


	// bind shaders
	m_pShaders->bind();
	BOOST_SCOPE_EXIT(m_pShaders) { m_pShaders->release(); } BOOST_SCOPE_EXIT_END
	LOGGLERR(pGl);

	if(m_bLightsNeedUpdate) UpdateLights();
	if(m_bBTrafoNeedsUpdate) UpdateBTrafo();

	// set cam matrix
	m_pShaders->setUniformValue(m_uniMatrixCam, m_matCam);
	m_pShaders->setUniformValue(m_uniMatrixCamInv, m_matCam_inv);


	auto colOverride = m::create<t_vec_gl>({1,1,1,1});
	auto colHighlight = m::create<t_vec_gl>({1,1,1,1});

	// render triangle geometry
	for(const auto& obj : m_objs)
	{
		const GlPlotObj *linkedObj = &obj;
		if(obj.linkedObj)
		{
			// get linked object
			linkedObj = &m_objs[*obj.linkedObj];

			// override constant color for linked object
			if(obj.m_highlighted)
				m_pShaders->setUniformValue(m_uniConstCol, colHighlight);
			else
				m_pShaders->setUniformValue(m_uniConstCol, obj.m_color);
		}
		else
		{
			// set override color to white for non-linked objects
			m_pShaders->setUniformValue(m_uniConstCol, colOverride);
		}

		if(!obj.m_visible || !obj.m_valid) continue;


		m_pShaders->setUniformValue(m_uniMatrixObj, obj.m_mat);

		// set to untransformed coordinate system if the object is invariant
		m_pShaders->setUniformValue(m_uniCoordSys, linkedObj->m_invariant ? 0 : m_iCoordSys.load());


		// main vertex array object
		pGl->glBindVertexArray(linkedObj->m_vertexarr);


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
void GlPlot_impl::DoPaintNonGL(QPainter &painter)
{
	const t_mat_gl matUnit = m::unit<t_mat_gl>();

	QFont fontOrig = painter.font();
	QPen penOrig = painter.pen();

	QPen penLabel(Qt::black);
	painter.setPen(penLabel);


	// coordinate labels
	painter.drawText(GlToScreenCoords(m::create<t_vec_gl>({0.,0.,0.,1.})), "0");
	for(t_real_gl f=-std::floor(m_CoordMax); f<=std::floor(m_CoordMax); f+=0.5)
	{
		if(m::equals<t_real_gl>(f, 0))
			continue;

		std::ostringstream ostrF;
		ostrF << f;
		painter.drawText(GlToScreenCoords(m::create<t_vec_gl>({f,0.,0.,1.})), ostrF.str().c_str());
		painter.drawText(GlToScreenCoords(m::create<t_vec_gl>({0.,f,0.,1.})), ostrF.str().c_str());
		painter.drawText(GlToScreenCoords(m::create<t_vec_gl>({0.,0.,f,1.})), ostrF.str().c_str());
	}

	painter.drawText(GlToScreenCoords(m::create<t_vec_gl>({m_CoordMax*t_real_gl(1.2), 0., 0., 1.})), "x");
	painter.drawText(GlToScreenCoords(m::create<t_vec_gl>({0., m_CoordMax*t_real_gl(1.2), 0., 1.})), "y");
	painter.drawText(GlToScreenCoords(m::create<t_vec_gl>({0., 0., m_CoordMax*t_real_gl(1.2), 1.})), "z");


	// render object labels
	for(const auto& obj : m_objs)
	{
		if(!obj.m_visible || !obj.m_valid) continue;

		if(obj.m_label != "")
		{
			const t_mat_gl *coordTrafo = &matUnit;
			if(m_iCoordSys == 1 && !obj.m_invariant) coordTrafo = &m_matA;

			t_vec3_gl posLabel3d = (*coordTrafo) * obj.m_mat * obj.m_labelPos;
			auto posLabel2d = GlToScreenCoords(m::create<t_vec_gl>({posLabel3d[0], posLabel3d[1], posLabel3d[2], 1.}));

			QFont fontLabel = fontOrig;
			QPen penLabel = penOrig;

			fontLabel.setStyleStrategy(QFont::StyleStrategy(QFont::OpenGLCompatible | QFont::PreferAntialias | QFont::PreferQuality));
			fontLabel.setWeight(QFont::Medium);
			//penLabel.setColor(QColor(int((1.-obj.m_color[0])*255.), int((1.-obj.m_color[1])*255.), int((1.-obj.m_color[2])*255.), int(obj.m_color[3]*255.)));
			penLabel.setColor(QColor(0,0,0,255));
			painter.setFont(fontLabel);
			painter.setPen(penLabel);
			painter.drawText(posLabel2d, obj.m_label.c_str());

			fontLabel.setWeight(QFont::Normal);
			penLabel.setColor(QColor(int(obj.m_color[0]*255.), int(obj.m_color[1]*255.), int(obj.m_color[2]*255.), int(obj.m_color[3]*255.)));
			painter.setFont(fontLabel);
			painter.setPen(penLabel);
			painter.drawText(posLabel2d, obj.m_label.c_str());
		}
	}

	// restore original styles
	painter.setFont(fontOrig);
	painter.setPen(penOrig);
}


void GlPlot_impl::paintGL()
{
	if(!m_bPlatformSupported) return;
	QMutexLocker _locker{&m_mutexObj};

	if constexpr(!m_isthreaded)
	{
		if(auto *pContext = m_pPlot->context(); !pContext) return;
		QPainter painter(m_pPlot);
		painter.setRenderHint(QPainter::HighQualityAntialiasing);

		// gl painting
		{
			BOOST_SCOPE_EXIT(m_pPlot, &painter) { painter.endNativePainting(); } BOOST_SCOPE_EXIT_END

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

		QMetaObject::invokeMethod(m_pPlot, &GlPlot::MoveContextToThread, Qt::ConnectionType::BlockingQueuedConnection);
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
				m_pPlot->GetImpl()->RequestPlotUpdate();
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

GlPlot::GlPlot(QWidget *pParent) : QOpenGLWidget(pParent),
	m_thread_impl(std::make_unique<QThread>(this)),
	m_impl(std::make_unique<GlPlot_impl>(this))
{
	qRegisterMetaType<std::size_t>("std::size_t");

	if constexpr(m_isthreaded)
	{
		m_impl->moveToThread(m_thread_impl.get());

		connect(m_thread_impl.get(), &QThread::started, m_impl.get(), &GlPlot_impl::startedThread);
		connect(m_thread_impl.get(), &QThread::finished, m_impl.get(), &GlPlot_impl::stoppedThread);
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
		m_impl->initialiseGL();
		if(m_impl->IsInitialised())
			emit AfterGLInitialisation();
	}
}


void GlPlot::resizeGL(int w, int h)
{
	if constexpr(!m_isthreaded)
	{
		m_impl->SetScreenDims(w, h);
		m_impl->resizeGL();
	}
}


void GlPlot::paintGL()
{
	if constexpr(!m_isthreaded)
	{
		m_impl->paintGL();
	}
}


void GlPlot::mouseMoveEvent(QMouseEvent *pEvt)
{
	m_impl->mouseMoveEvent(pEvt->localPos());
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
		m_impl->ResetZoom();
	if(m_mouseDown[2])
		m_impl->BeginRotation();

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
		m_impl->EndRotation();

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
	m_impl->zoom(degrees);

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
		QMetaObject::invokeMethod(m_impl.get(), &GlPlot_impl::paintGL, Qt::ConnectionType::QueuedConnection);
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
		m_impl->SetScreenDims(w, h);
	}
}

// ----------------------------------------------------------------------------
