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

#ifndef __MAG_GL_COMMON_H__
#define __MAG_GL_COMMON_H__

#if QT_VERSION >= 0x060000
	#include <QtOpenGL/QOpenGLShaderProgram>
	#include <QtOpenGL/QOpenGLVertexArrayObject>
	#include <QtOpenGL/QOpenGLBuffer>
	#include <QtOpenGL/QOpenGLFramebufferObject>
	#include <QtOpenGL/QOpenGLFramebufferObjectFormat>
	#include <QtOpenGL/QOpenGLPaintDevice>
	#include <QtOpenGLWidgets/QOpenGLWidget>

	#include <QtCore/QRecursiveMutex>
	using t_qt_mutex = QRecursiveMutex;
#else
	#include <QtGui/QOpenGLShaderProgram>
	#include <QtGui/QOpenGLVertexArrayObject>
	#include <QtGui/QOpenGLBuffer>
	#include <QtGui/QOpenGLFramebufferObject>
	#include <QtGui/QOpenGLFramebufferObjectFormat>
	#include <QtGui/QOpenGLPaintDevice>
	#include <QtWidgets/QOpenGLWidget>

	#include <QtCore/QMutex>
	using t_qt_mutex = QMutex;
#endif

#include <QtGui/QMatrix4x4>
#include <QtGui/QVector4D>
#include <QtGui/QVector3D>
#include <QtGui/QVector2D>

#include <memory>

#include "../maths.h"



// ----------------------------------------------------------------------------
// GL versions
#if !defined(_GL_MAJ_VER) || !defined(_GL_MIN_VER)
	#define _GL_MAJ_VER 3
	#define _GL_MIN_VER 3
#endif

#if _GL_MAJ_VER<=3 && _GL_MIN_VER<2
	#if !defined(_GL_SUFFIX)
		#define _GL_SUFFIX
	#endif

	#if _GL_MAJ_VER==3 && _GL_MIN_VER==1
		#define _GLSL_MAJ_VER 1
		#define _GLSL_MIN_VER 4
	#elif _GL_MAJ_VER==3 && _GL_MIN_VER==0
		#define _GLSL_MAJ_VER 1
		#define _GLSL_MIN_VER 3
	#elif _GL_MAJ_VER==2 && _GL_MIN_VER==1
		#define _GLSL_MAJ_VER 1
		#define _GLSL_MIN_VER 2
	#elif _GL_MAJ_VER==2 && _GL_MIN_VER==0
		#define _GLSL_MAJ_VER 1
		#define _GLSL_MIN_VER 1
	#endif
#else
	#if !defined(_GL_SUFFIX)
		#define _GL_SUFFIX _Core
	#endif

	#if _GL_MAJ_VER==3 && _GL_MIN_VER==2
		#define _GLSL_MAJ_VER 1
		#define _GLSL_MIN_VER 5
	#else
		#define _GLSL_MAJ_VER _GL_MAJ_VER
		#define _GLSL_MIN_VER _GL_MIN_VER
	#endif
#endif


// GL functions include
#if QT_VERSION >= 0x060000
	#define _GL_INC_IMPL(MAJ, MIN, SUFF) <QtOpenGL/QOpenGLFunctions_ ## MAJ ## _ ## MIN ## SUFF>
#else
	#define _GL_INC_IMPL(MAJ, MIN, SUFF) <QtGui/QOpenGLFunctions_ ## MAJ ## _ ## MIN ## SUFF>
#endif
#define _GL_INC(MAJ, MIN, SUFF) _GL_INC_IMPL(MAJ, MIN, SUFF)
#include _GL_INC(_GL_MAJ_VER, _GL_MIN_VER, _GL_SUFFIX)

// GL functions typedef
#define _GL_FUNC_IMPL(MAJ, MIN, SUFF) QOpenGLFunctions_ ## MAJ ## _ ## MIN ## SUFF
#define _GL_FUNC(MAJ, MIN, SUFF) _GL_FUNC_IMPL(MAJ, MIN, SUFF)
using qgl_funcs = _GL_FUNC(_GL_MAJ_VER, _GL_MIN_VER, _GL_SUFFIX);


// GL error codes: https://www.khronos.org/opengl/wiki/OpenGL_Error
#define LOGGLERR(pGl) { while(true) {	\
		auto err = pGl->glGetError();	\
		if(err == GL_NO_ERROR) break;	\
		std::cerr << "GL error in " << __func__ \
			<< ", file: " << __FILE__ \
			<< ", line " << std::dec <<  __LINE__  \
			<< ": " << std::hex << "0x" << err \
			<< "." << std::endl; \
	}}
// ----------------------------------------------------------------------------


namespace tl2 {

// ----------------------------------------------------------------------------
// types
using t_real_gl = GLfloat;
//using t_real_gl = GLdouble;
using t_vec2_gl = tl2::qvecN_adapter<int, 2, t_real_gl, QVector2D>;
using t_vec3_gl = tl2::qvecN_adapter<int, 3, t_real_gl, QVector3D>;
using t_vec_gl = tl2::qvecN_adapter<int, 4, t_real_gl, QVector4D>;
using t_mat_gl = tl2::qmatNN_adapter<int, 4, 4, t_real_gl, QMatrix4x4>;
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// plotter objects
enum class GlRenderObjType
{
	TRIANGLES,
	LINES
};


#ifndef GlPlotObjType
	#define GlPlotObjType GlRenderObjType
#endif


struct GlRenderObj
{
	GlPlotObjType m_type = GlPlotObjType::TRIANGLES;

	std::shared_ptr<QOpenGLVertexArrayObject> m_pvertexarr{};
	std::shared_ptr<QOpenGLBuffer> m_pvertexbuf{};
	std::shared_ptr<QOpenGLBuffer> m_pnormalsbuf{};
	std::shared_ptr<QOpenGLBuffer> m_puvbuf{};
	std::shared_ptr<QOpenGLBuffer> m_pcolourbuf{};

	std::vector<t_vec3_gl> m_vertices{}, m_triangles{}, m_uvs{};

	t_vec_gl m_colour = tl2::create<t_vec_gl>({ 0., 0., 1., 1. });	// rgba
};


struct GlPlotObj : public GlRenderObj
{
	// does not define a geometry itself, but just links to another object
	std::optional<std::size_t> linkedObj{};

	t_mat_gl m_mat = tl2::unit<t_mat_gl>();

	bool m_invariant = false;	// invariant to A, B matrices
	bool m_visible = true;		// object shown?
	bool m_highlighted = false;	// object highlighted?
	bool m_valid = true;		// object deleted?

	t_vec3_gl m_labelPos = tl2::create<t_vec3_gl>({0., 0., 0.});
	std::string m_label{};
	std::string m_datastr{};

	t_vec3_gl m_boundingSpherePos = tl2::create<t_vec3_gl>({ 0., 0., 0. });
	t_real_gl m_boundingSphereRad = 0.;
};
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// functions
// GL surface format
extern void set_gl_format(bool bCore = true, int iMajorVer = 3, int iMinorVer = 3, int iSamples = 8);
extern QSurfaceFormat gl_format(
	bool bCore = true, int iMajorVer = 3, int iMinorVer = 3, int iSamples = 8,
	QSurfaceFormat surf = QSurfaceFormat::defaultFormat());

// get gl functions
extern qgl_funcs* get_gl_functions(QOpenGLWidget *pGLWidget);

// create a triangle-based object
extern bool create_triangle_object(QOpenGLWidget* pGLWidget, GlRenderObj& obj,
	const std::vector<t_vec3_gl>& verts, const std::vector<t_vec3_gl>& triagverts,
	const std::vector<t_vec3_gl>& norms, const std::vector<t_vec3_gl>& uvs, const t_vec_gl& colour,
	bool bUseVertsAsNorm, GLint attrVertex, GLint attrVertexNormal, GLint attrVertexcolour,
	GLint attrTextureCoords=-1);

// create a line-based object
extern bool create_line_object(QOpenGLWidget* pGLWidget, GlRenderObj& obj,
	const std::vector<t_vec3_gl>& verts, const t_vec_gl& colour,
	GLint attrVertex, GLint attrVertexcolour);


extern void delete_render_object(GlRenderObj& obj);
// ----------------------------------------------------------------------------

}
#endif
