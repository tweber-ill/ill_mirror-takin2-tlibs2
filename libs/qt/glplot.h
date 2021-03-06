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
 */

#ifndef __MAG_GL_PLOT_H__
#define __MAG_GL_PLOT_H__

#include <QtCore/QThread>
#include <QtCore/QTimer>
#include <QtWidgets/QDialog>
#include <QtGui/QMouseEvent>

#include <memory>
#include <chrono>
#include <atomic>

#include "gl.h"
#include "../maths.h"


namespace tl2 {

// ----------------------------------------------------------------------------
class GlPlot;


/**
 * GL plot renderer
 */
class GlPlotRenderer : public QObject
{ Q_OBJECT
private:
#if QT_VERSION >= 0x060000
	t_qt_mutex m_mutexObj;
#else
	t_qt_mutex m_mutexObj{QMutex::Recursive};
#endif


protected:
	GlPlot *m_pPlot = nullptr;
	std::string m_strGlVer, m_strGlShaderVer, m_strGlVendor, m_strGlRenderer;

	std::shared_ptr<QOpenGLShaderProgram> m_pShaders;

	GLint m_attrVertex = -1;
	GLint m_attrVertexNorm = -1;
	GLint m_attrVertexCol = -1;
	GLint m_uniConstCol = -1;
	GLint m_uniLightPos = -1;
	GLint m_uniNumActiveLights = -1;
	GLint m_uniMatrixProj = -1;
	GLint m_uniMatrixCam = -1;
	GLint m_uniMatrixCamInv = -1;
	GLint m_uniMatrixObj = -1;
	GLint m_uniMatrixA = -1;
	GLint m_uniMatrixB = -1;
	GLint m_uniCoordSys = -1;

	t_mat_gl m_matPerspective = tl2::unit<t_mat_gl>();
	t_mat_gl m_matPerspective_inv = tl2::unit<t_mat_gl>();
	t_mat_gl m_matViewport = tl2::unit<t_mat_gl>();
	t_mat_gl m_matViewport_inv = tl2::unit<t_mat_gl>();
	t_mat_gl m_matCamBase = tl2::create<t_mat_gl>({1,0,0,0,  0,1,0,0,  0,0,1,-5,  0,0,0,1});
	t_mat_gl m_matCamRot = tl2::unit<t_mat_gl>();
	t_mat_gl m_matCam = tl2::unit<t_mat_gl>();
	t_mat_gl m_matCam_inv = tl2::unit<t_mat_gl>();
	t_mat_gl m_matA = tl2::unit<t_mat_gl>();
	t_mat_gl m_matB = tl2::unit<t_mat_gl>();
	t_vec_gl m_vecCamX = tl2::create<t_vec_gl>({1.,0.,0.,0.});
	t_vec_gl m_vecCamY = tl2::create<t_vec_gl>({0.,1.,0.,0.});
	t_real_gl m_phi_saved = 0, m_theta_saved = 0;
	t_real_gl m_zoom = 1.;
	t_real_gl m_CoordMax = 2.5;		// extent of coordinate axes

	std::atomic<bool> m_bPlatformSupported = true;
	std::atomic<bool> m_bInitialised = false;
	std::atomic<bool> m_bWantsResize = false;
	std::atomic<bool> m_bPickerEnabled = true;
	std::atomic<bool> m_bPickerNeedsUpdate = false;
	std::atomic<bool> m_bLightsNeedUpdate = false;
	std::atomic<bool> m_bBTrafoNeedsUpdate = false;
	std::atomic<int> m_iCoordSys = 0;
	std::atomic<int> m_iScreenDims[2] = { 800, 600 };
	t_real_gl m_pickerSphereRadius = 1;

	std::vector<t_vec3_gl> m_lights;
	std::vector<GlPlotObj> m_objs;

	QPointF m_posMouse;
	QPointF m_posMouseRotationStart, m_posMouseRotationEnd;
	bool m_bInRotation = false;

	QTimer m_timer;


protected:
	qgl_funcs* GetGlFunctions() { return get_gl_functions((QOpenGLWidget*)m_pPlot); }

	void UpdateCam();
	void RequestPlotUpdate();
	void UpdatePicker();
	void UpdateLights();
	void UpdateBTrafo();

	void DoPaintGL(qgl_funcs *pGL);
	void DoPaintNonGL(QPainter &painter);

	void tick(const std::chrono::milliseconds& ms);


public:
	GlPlotRenderer(GlPlot *pPlot = nullptr);
	virtual ~GlPlotRenderer();

	static constexpr bool m_isthreaded = false;
	static constexpr bool m_usetimer = false;

	std::tuple<std::string, std::string, std::string, std::string>
		GetGlDescr() const { return std::make_tuple(m_strGlVer, m_strGlShaderVer, m_strGlVendor, m_strGlRenderer); }

	QPointF GlToScreenCoords(const t_vec_gl& vec, bool *pVisible=nullptr) const;

	void SetCamBase(const t_mat_gl& mat, const t_vec_gl& vecX, const t_vec_gl& vecY)
	{ m_matCamBase = mat; m_vecCamX = vecX; m_vecCamY = vecY; UpdateCam(); }
	void SetPickerSphereRadius(t_real_gl rad) { m_pickerSphereRadius = rad; }

	GlPlotObj CreateTriangleObject(const std::vector<t_vec3_gl>& verts,
		const std::vector<t_vec3_gl>& triag_verts, const std::vector<t_vec3_gl>& norms,
		const t_vec_gl& colour, bool bUseVertsAsNorm=false);
	GlPlotObj CreateLineObject(const std::vector<t_vec3_gl>& verts, const t_vec_gl& colour);

	std::size_t GetNumObjects() const { return m_objs.size(); }
	void RemoveObject(std::size_t obj);
	std::size_t AddLinkedObject(std::size_t linkTo,
		t_real_gl x=0, t_real_gl y=0, t_real_gl z=0,
		t_real_gl r=1, t_real_gl g=1, t_real_gl b=1, t_real_gl a=1);
	std::size_t AddSphere(t_real_gl rad=1,
		t_real_gl x=0, t_real_gl y=0, t_real_gl z=0,
		t_real_gl r=0, t_real_gl g=0, t_real_gl b=0, t_real_gl a=1);
	std::size_t AddCylinder(t_real_gl rad=1, t_real_gl h=1,
		t_real_gl x=0, t_real_gl y=0, t_real_gl z=0,
		t_real_gl r=0, t_real_gl g=0, t_real_gl b=0, t_real_gl a=1);
	std::size_t AddCone(t_real_gl rad=1, t_real_gl h=1,
		t_real_gl x=0, t_real_gl y=0, t_real_gl z=0,
		t_real_gl r=0, t_real_gl g=0, t_real_gl b=0, t_real_gl a=1);
	std::size_t AddArrow(t_real_gl rad=1, t_real_gl h=1,
		t_real_gl x=0, t_real_gl y=0, t_real_gl z=0,
		t_real_gl r=0, t_real_gl g=0, t_real_gl b=0, t_real_gl a=1);
	std::size_t AddTriangleObject(const std::vector<t_vec3_gl>& triag_verts,
		const std::vector<t_vec3_gl>& triag_norms,
		t_real_gl r=0, t_real_gl g=0, t_real_gl b=0, t_real_gl a=1);
	std::size_t AddCoordinateCross(t_real_gl min, t_real_gl max);

	void SetObjectMatrix(std::size_t idx, const t_mat_gl& mat);
	void SetObjectCol(std::size_t idx, t_real_gl r, t_real_gl g, t_real_gl b, t_real_gl a=1);
	void SetObjectLabel(std::size_t idx, const std::string& label);
	void SetObjectDataString(std::size_t idx, const std::string& data);
	void SetObjectVisible(std::size_t idx, bool visible);
	void SetObjectHighlight(std::size_t idx, bool highlight);

	const std::string& GetObjectLabel(std::size_t idx) const;
	const std::string& GetObjectDataString(std::size_t idx) const;
	bool GetObjectHighlight(std::size_t idx) const;

	void SetScreenDims(int w, int h);
	void SetCoordMax(t_real_gl d) { m_CoordMax = d; }

	void SetLight(std::size_t idx, const t_vec3_gl& pos);

	void SetBTrafo(const t_mat_gl& matB, const t_mat_gl* matA = nullptr);
	void SetCoordSys(int iSys);

	bool IsInitialised() const { return m_bInitialised; }


public slots:
	void paintGL();

	void startedThread();
	void stoppedThread();

	void initialiseGL();
	void resizeGL();

	void mouseMoveEvent(const QPointF& pos);
	void zoom(t_real_gl val);
	void ResetZoom();

	void BeginRotation();
	void EndRotation();

	void EnablePicker(bool b);


protected slots:
	void tick();

signals:
	void PickerIntersection(const t_vec3_gl* pos, std::size_t objIdx, const t_vec3_gl* posSphere);
};



/**
 * GL plotter widget
 */
class GlPlot : public QOpenGLWidget
{ Q_OBJECT
public:
	GlPlot(QWidget *pParent = nullptr);
	virtual ~GlPlot();

	GlPlotRenderer* GetRenderer() { return m_renderer.get(); }
	static constexpr bool m_isthreaded = GlPlotRenderer::m_isthreaded;

protected:
	virtual void paintEvent(QPaintEvent*) override;
	virtual void initializeGL() override;
	virtual void paintGL() override;
	virtual void resizeGL(int w, int h) override;

	virtual void mouseMoveEvent(QMouseEvent *pEvt) override;
	virtual void mousePressEvent(QMouseEvent *Evt) override;
	virtual void mouseReleaseEvent(QMouseEvent *Evt) override;
	virtual void wheelEvent(QWheelEvent *pEvt) override;

private:
#if QT_VERSION >= 0x060000
	mutable t_qt_mutex m_mutex;
#else
	mutable t_qt_mutex m_mutex{QMutex::Recursive};
#endif

	std::unique_ptr<GlPlotRenderer> m_renderer;
	std::unique_ptr<QThread> m_thread_impl;
	bool m_mouseMovedBetweenDownAndUp = 0;
	bool m_mouseDown[3] = {0,0,0};

public:
	t_qt_mutex* GetMutex() { return &m_mutex; }

	void MoveContextToThread();
	bool IsContextInThread() const;

protected slots:
	void beforeComposing();
	void afterComposing();
	void beforeResizing();
	void afterResizing();

signals:
	void AfterGLInitialisation();
	void GLInitialisationFailed();

	void MouseDown(bool left, bool mid, bool right);
	void MouseUp(bool left, bool mid, bool right);
	void MouseClick(bool left, bool mid, bool right);
};
// ----------------------------------------------------------------------------

}
#endif
