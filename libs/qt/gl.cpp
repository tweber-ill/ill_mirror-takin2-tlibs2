/**
 * tlibs2 -- common gl functions
 * @author Tobias Weber <tweber@ill.fr>
 * @date 2017-2021
 * @note The present version was forked on 8-Nov-2018 from my privately developed "magtools" project (https://github.com/t-weber/magtools).
 * @license GPLv3, see 'LICENSE' file
 *
 * References:
 *   - http://doc.qt.io/qt-5/qopenglwidget.html#details
 *   - http://code.qt.io/cgit/qt/qtbase.git/tree/examples/opengl/threadedqopenglwidget
 */

#include "gl.h"

#include <QtGui/QOpenGLContext>
#include <QtGui/QOpenGLFunctions>
#include <QtGui/QSurfaceFormat>

#if QT_VERSION >= 0x060000
	#include <QtOpenGL/QOpenGLVersionFunctionsFactory>
#endif

#include <iostream>
#include <boost/scope_exit.hpp>
#include <boost/preprocessor/stringize.hpp>


#pragma message("Compiling for GL version " BOOST_PP_STRINGIZE(_GL_MAJ_VER) "." BOOST_PP_STRINGIZE(_GL_MIN_VER) \
	" and GLSL version " BOOST_PP_STRINGIZE(_GLSL_MAJ_VER) BOOST_PP_STRINGIZE(_GLSL_MIN_VER) "0.")


namespace tl2 {
// ----------------------------------------------------------------------------
// functions
// ----------------------------------------------------------------------------

/**
 * create a gl surface format
 */
extern QSurfaceFormat gl_format(
	bool bCore, int iMajorVer, int iMinorVer, int iSamples, QSurfaceFormat surf)
{
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

	return surf;
}


/**
 * set the default gl surface format
 */
void set_gl_format(bool bCore, int iMajorVer, int iMinorVer, int iSamples)
{
	QSurfaceFormat surf = gl_format(bCore, iMajorVer, iMinorVer, iSamples, QSurfaceFormat::defaultFormat());
	QSurfaceFormat::setDefaultFormat(surf);
}


/**
 * return gl functions for current version
 */
qgl_funcs* get_gl_functions(QOpenGLWidget *pGLWidget)
{
	qgl_funcs *pGl = nullptr;

	if constexpr(std::is_same_v<qgl_funcs, QOpenGLFunctions>)
	{
		pGl = (qgl_funcs*)pGLWidget->context()->functions();
	}
	else
	{
#if QT_VERSION >= 0x060000
		pGl = QOpenGLVersionFunctionsFactory::get<qgl_funcs>(pGLWidget->context());
#else
		pGl = (qgl_funcs*)pGLWidget->context()->versionFunctions<qgl_funcs>();
#endif
	}

	if(!pGl)
		std::cerr << "No suitable GL interface found." << std::endl;

	return pGl;
}


/**
 * creates a triangle-based 3d object
 */
bool create_triangle_object(QOpenGLWidget* pGLWidget, GlRenderObj& obj,
	const std::vector<t_vec3_gl>& verts, const std::vector<t_vec3_gl>& triagverts,
	const std::vector<t_vec3_gl>& norms, const std::vector<t_vec3_gl>& uvs,
	const t_vec_gl& colour, bool bUseVertsAsNorm,
	GLint attrVertex, GLint attrVertexNormal, GLint attrVertexcolour, GLint attrTextureCoords)
{
	// TODO: move context to calling thread
	pGLWidget->makeCurrent();
	BOOST_SCOPE_EXIT(pGLWidget) { pGLWidget->doneCurrent(); } BOOST_SCOPE_EXIT_END

	qgl_funcs* pGl = get_gl_functions(pGLWidget);
	if(!pGl) return false;

	obj.m_type = GlRenderObjType::TRIANGLES;
	obj.m_colour = colour;

	// flatten vertex array into raw float array
	auto to_float_array = [](const std::vector<t_vec3_gl>& verts, 
		int iRepeat=1, int iElems=3, bool bNorm=false, t_real_gl lastElem=1.)
		-> std::vector<t_real_gl>
	{
		std::vector<t_real_gl> vecRet;
		vecRet.reserve(iRepeat*verts.size()*iElems);

		for(const t_vec3_gl& vert : verts)
		{
			t_real_gl norm = bNorm ? tl2::norm<t_vec3_gl>(vert) : 1;

			for(int i=0; i<iRepeat; ++i)
			{
				for(int iElem=0; iElem<iElems; ++iElem)
				{
					if(iElem < vert.size())
						vecRet.push_back(vert[iElem] / norm);
					else
						vecRet.push_back(lastElem);
				}
			}
		}

		return vecRet;
	};

	// main vertex array object
	pGl->glGenVertexArrays(1, &obj.m_vertexarr);
	pGl->glBindVertexArray(obj.m_vertexarr);

	// vertices
	if(attrVertex >= 0)
	{
		obj.m_pvertexbuf = std::make_shared<QOpenGLBuffer>(QOpenGLBuffer::VertexBuffer);

		if(!obj.m_pvertexbuf->create())
			std::cerr << "Cannot create vertex buffer." << std::endl;
		if(!obj.m_pvertexbuf->bind())
			std::cerr << "Cannot bind vertex buffer." << std::endl;
		BOOST_SCOPE_EXIT(&obj) { obj.m_pvertexbuf->release(); } BOOST_SCOPE_EXIT_END

		auto vecVerts = to_float_array(triagverts, 1, 4, false, 1.);
		obj.m_pvertexbuf->allocate(vecVerts.data(), vecVerts.size()*sizeof(typename decltype(vecVerts)::value_type));
		pGl->glVertexAttribPointer(attrVertex, 4, GL_FLOAT, 0, 0, nullptr);
	}

	// normals
	if(attrVertexNormal >= 0)
	{
		obj.m_pnormalsbuf = std::make_shared<QOpenGLBuffer>(QOpenGLBuffer::VertexBuffer);

		obj.m_pnormalsbuf->create();
		obj.m_pnormalsbuf->bind();
		BOOST_SCOPE_EXIT(&obj) { obj.m_pnormalsbuf->release(); } BOOST_SCOPE_EXIT_END

		auto vecNorms = bUseVertsAsNorm ? to_float_array(triagverts, 1, 4, true, 0.) : to_float_array(norms, 3, 4, false, 0.);
		obj.m_pnormalsbuf->allocate(vecNorms.data(), vecNorms.size()*sizeof(typename decltype(vecNorms)::value_type));
		pGl->glVertexAttribPointer(attrVertexNormal, 4, GL_FLOAT, 0, 0, nullptr);
	}

	// colours
	if(attrVertexcolour >= 0)
	{
		obj.m_pcolourbuf = std::make_shared<QOpenGLBuffer>(QOpenGLBuffer::VertexBuffer);

		obj.m_pcolourbuf->create();
		obj.m_pcolourbuf->bind();
		BOOST_SCOPE_EXIT(&obj) { obj.m_pcolourbuf->release(); } BOOST_SCOPE_EXIT_END

		std::vector<t_real_gl> vecCols;
		vecCols.reserve(4*triagverts.size());
		for(std::size_t iVert=0; iVert<triagverts.size(); ++iVert)
		{
			for(int icol=0; icol<obj.m_colour.size(); ++icol)
				vecCols.push_back(obj.m_colour[icol]);
		}

		obj.m_pcolourbuf->allocate(vecCols.data(), vecCols.size()*sizeof(typename decltype(vecCols)::value_type));
		pGl->glVertexAttribPointer(attrVertexcolour, 4, GL_FLOAT, 0, 0, nullptr);
	}

	// texture uv coordinates
	if(attrTextureCoords >= 0)
	{
		obj.m_puvbuf = std::make_shared<QOpenGLBuffer>(QOpenGLBuffer::VertexBuffer);

		obj.m_puvbuf->create();
		obj.m_puvbuf->bind();
		BOOST_SCOPE_EXIT(&obj) { obj.m_puvbuf->release(); } BOOST_SCOPE_EXIT_END

		auto vecUVs = to_float_array(uvs, 1, 2);
		obj.m_puvbuf->allocate(vecUVs.data(), vecUVs.size()*sizeof(typename decltype(vecUVs)::value_type));
		pGl->glVertexAttribPointer(attrTextureCoords, 2, GL_FLOAT, 0, 0, nullptr);
	}


	obj.m_vertices = std::move(verts);
	obj.m_triangles = std::move(triagverts);
	obj.m_uvs = std::move(uvs);
	LOGGLERR(pGl)

	return true;
}


/**
 * creates a line-based 3d object
 */
bool create_line_object(QOpenGLWidget* pGLWidget, GlRenderObj& obj,
	const std::vector<t_vec3_gl>& verts, const t_vec_gl& colour,
	GLint attrVertex, GLint attrVertexcolour)
{
	// TODO: move context to calling thread
	pGLWidget->makeCurrent();
	BOOST_SCOPE_EXIT(pGLWidget) { pGLWidget->doneCurrent(); } BOOST_SCOPE_EXIT_END

	qgl_funcs* pGl = get_gl_functions(pGLWidget);
	if(!pGl) return false;

	//GLint attrVertex = m_attrVertex;
	//GLint attrVertexcolour = m_attrVertexCol;

	obj.m_type = GlRenderObjType::LINES;
	obj.m_colour = colour;

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

	{	// colours
		obj.m_pcolourbuf = std::make_shared<QOpenGLBuffer>(QOpenGLBuffer::VertexBuffer);

		obj.m_pcolourbuf->create();
		obj.m_pcolourbuf->bind();
		BOOST_SCOPE_EXIT(&obj) { obj.m_pcolourbuf->release(); } BOOST_SCOPE_EXIT_END

		std::vector<t_real_gl> vecCols;
		vecCols.reserve(4*verts.size());
		for(std::size_t iVert=0; iVert<verts.size(); ++iVert)
		{
			for(int icol=0; icol<obj.m_colour.size(); ++icol)
				vecCols.push_back(obj.m_colour[icol]);
		}

		obj.m_pcolourbuf->allocate(vecCols.data(), vecCols.size()*sizeof(typename decltype(vecCols)::value_type));
		pGl->glVertexAttribPointer(attrVertexcolour, 4, GL_FLOAT, 0, 0, nullptr);
	}


	obj.m_vertices = std::move(verts);
	LOGGLERR(pGl)

	return true;
}
// ----------------------------------------------------------------------------

}
