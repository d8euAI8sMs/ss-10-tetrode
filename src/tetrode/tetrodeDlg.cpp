// tetrodeDlg.cpp : implementation file
//

#include "stdafx.h"
#include "tetrode.h"
#include "tetrodeDlg.h"
#include "afxdialogex.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// CTetrodeDlg dialog



CTetrodeDlg::CTetrodeDlg(CWnd* pParent /*=NULL*/)
    : CSimulationDialog(CTetrodeDlg::IDD, pParent)
    , data(model::make_model_data())
    , m_bPointsVisible(TRUE)
    , m_bTriangulationVisible(FALSE)
    , m_bDirichletCellsVisible(FALSE)
    , m_bIsolineVisible(FALSE)
    , m_nIsolineCount(100)
    , m_lfIsolineDelta(0.1)
    , m_lfAccuracyGoal(1e-4)
    , m_lfAnodeBeginPotential(-10)
    , m_lfAnodeEndPotential(10)
    , m_nAnodePotentialSamples(21)
{
    m_hIcon = AfxGetApp()->LoadIcon(IDR_MAINFRAME);
}

void CTetrodeDlg::DoDataExchange(CDataExchange* pDX)
{
    CSimulationDialog::DoDataExchange(pDX);
    DDX_Control(pDX, IDC_MODEL_PLOT, m_cSystemPlot);
    DDX_Control(pDX, IDC_ANODE_CURRENT_PLOT, m_cAnodeCurrentPlot);
    DDX_Text(pDX, IDC_EDIT1, data.params->v0);
    DDX_Text(pDX, IDC_EDIT2, data.params->i0);
    DDX_Text(pDX, IDC_EDIT2, data.params->i0);
    DDX_Text(pDX, IDC_EDIT7, data.params->u1);
    DDX_Text(pDX, IDC_EDIT8, data.params->u2);
    DDX_Text(pDX, IDC_EDIT9, data.params->ua);
    DDX_Text(pDX, IDC_EDIT10, data.params->x1);
    DDX_Text(pDX, IDC_EDIT11, data.params->x2);
    DDX_Text(pDX, IDC_EDIT12, data.params->w);
    DDX_Text(pDX, IDC_EDIT3, data.params->dt);
    DDX_Text(pDX, IDC_EDIT4, data.params->ndt);
    DDX_Text(pDX, IDC_EDIT5, data.params->dx);
    DDX_Text(pDX, IDC_EDIT6, data.params->dy);
    DDX_Check(pDX, IDC_CHECK1, m_bPointsVisible);
    DDX_Check(pDX, IDC_CHECK2, m_bTriangulationVisible);
    DDX_Check(pDX, IDC_CHECK3, m_bDirichletCellsVisible);
    DDX_Check(pDX, IDC_CHECK4, m_bIsolineVisible);
    DDX_Text(pDX, IDC_EDIT17, m_nIsolineCount);
    DDX_Text(pDX, IDC_EDIT18, m_lfIsolineDelta);
    DDX_Text(pDX, IDC_EDIT15, m_lfAccuracyGoal);
    DDX_Text(pDX, IDC_EDIT16, m_nAnodePotentialSamples);
    DDX_Text(pDX, IDC_EDIT13, m_lfAnodeBeginPotential);
    DDX_Text(pDX, IDC_EDIT14, m_lfAnodeEndPotential);
}

BEGIN_MESSAGE_MAP(CTetrodeDlg, CSimulationDialog)
    ON_WM_PAINT()
    ON_WM_QUERYDRAGICON()
    ON_BN_CLICKED(IDC_BUTTON1, &CTetrodeDlg::OnBnClickedButton1)
    ON_BN_CLICKED(IDC_BUTTON3, &CTetrodeDlg::OnBnClickedButton3)
    ON_BN_CLICKED(IDC_BUTTON2, &CTetrodeDlg::OnBnClickedButton2)
    ON_BN_CLICKED(IDC_CHECK1, &CTetrodeDlg::OnBnClickedCheck1)
    ON_BN_CLICKED(IDC_CHECK2, &CTetrodeDlg::OnBnClickedCheck1)
    ON_BN_CLICKED(IDC_CHECK3, &CTetrodeDlg::OnBnClickedCheck1)
    ON_BN_CLICKED(IDC_CHECK4, &CTetrodeDlg::OnBnClickedCheck1)
END_MESSAGE_MAP()


// CTetrodeDlg message handlers

BOOL CTetrodeDlg::OnInitDialog()
{
    CSimulationDialog::OnInitDialog();

    // Set the icon for this dialog.  The framework does this automatically
    //  when the application's main window is not a dialog
    SetIcon(m_hIcon, TRUE);            // Set big icon
    SetIcon(m_hIcon, FALSE);        // Set small icon

    // TODO: Add extra initialization here

    OnBnClickedCheck1();

    m_cSystemPlot.plot_layer.with(model::make_root_drawable(
        plot::make_world_mapper(data.system_data.world), {
        data.system_data.dirichlet_cell_plot,
        data.system_data.triangulation_plot,
        data.system_data.system_plot,
        data.isoline_data.plot,
        data.system_data.point_plot
    }));

    m_cAnodeCurrentPlot.plot_layer.with(model::make_root_drawable(
        plot::make_world_mapper(data.func_data.world), {
        data.func_data.plot
    }));

    m_cSystemPlot.background = plot::palette::brush();
    m_cSystemPlot.triple_buffered = true;
    m_cAnodeCurrentPlot.background = plot::palette::brush();
    m_cAnodeCurrentPlot.triple_buffered = true;

    m_cSystemPlot.RedrawBuffer();
    m_cSystemPlot.SwapBuffers();
    m_cSystemPlot.RedrawWindow();
    m_cAnodeCurrentPlot.RedrawBuffer();
    m_cAnodeCurrentPlot.SwapBuffers();
    m_cAnodeCurrentPlot.RedrawWindow();

    return TRUE;  // return TRUE  unless you set the focus to a control
}

// If you add a minimize button to your dialog, you will need the code below
//  to draw the icon.  For MFC applications using the document/view model,
//  this is automatically done for you by the framework.

void CTetrodeDlg::OnPaint()
{
    if (IsIconic())
    {
        CPaintDC dc(this); // device context for painting

        SendMessage(WM_ICONERASEBKGND, reinterpret_cast<WPARAM>(dc.GetSafeHdc()), 0);

        // Center icon in client rectangle
        int cxIcon = GetSystemMetrics(SM_CXICON);
        int cyIcon = GetSystemMetrics(SM_CYICON);
        CRect rect;
        GetClientRect(&rect);
        int x = (rect.Width() - cxIcon + 1) / 2;
        int y = (rect.Height() - cyIcon + 1) / 2;

        // Draw the icon
        dc.DrawIcon(x, y, m_hIcon);
    }
    else
    {
        CSimulationDialog::OnPaint();
    }
}

// The system calls this function to obtain the cursor to display while the user drags
//  the minimized window.
HCURSOR CTetrodeDlg::OnQueryDragIcon()
{
    return static_cast<HCURSOR>(m_hIcon);
}


void CTetrodeDlg::OnBnClickedButton1()
{
    m_eSimulationMode = sm_it;
    UpdateData(TRUE);
    StartSimulationThread();
}


void CTetrodeDlg::OnBnClickedButton3()
{
    m_eSimulationMode = sm_iv;
    UpdateData(TRUE);
    StartSimulationThread();
}


void CTetrodeDlg::OnBnClickedButton2()
{
    StopSimulationThread();
}


void CTetrodeDlg::OnSimulation()
{
    omp_set_num_threads(omp_get_num_procs());
    simulation_mode sm;
    Invoke([&] () { sm = m_eSimulationMode; });
    model::adjust(*data.params, *data.system_data.world);
    model::update_system_data(*data.params, data.system_data);
    model::particle_particle pp
    (
        *data.params,
        data.system_data.mesh,
        *data.system_data.data
    );
    pp.init();
    model::finel_galerkin g
    (
        *data.params,
        data.system_data.mesh,
        pp.charges, m_lfAccuracyGoal, 1000
    );
    g.init();

    size_t ndt = 0;
    double current = 0;

    while (m_bWorking)
    {
        if (sm == sm_it)
        {
            g.update();

            g.next(pp.potential);

            pp.generate_particles();

            double q0 = pp.collect_charges_and_adjust_particles();

            if (m_bIsolineVisible)
            {
                model::find_isolines(*data.system_data.mesh,
                                     pp.potential,
                                     m_lfIsolineDelta,
                                     m_nIsolineCount,
                                     *data.isoline_data.data);
            }

            current += q0;
            ++ndt;

            if ((ndt % data.params->ndt) == 0)
            {
                current /= ndt * data.params->dt;
                double t = 0;
                if (!data.func_data.data->empty())
                    t = data.func_data.data->back().x + data.params->dt;
                data.func_data.world->xmax = t;
                data.func_data.world->ymax =
                    max(data.func_data.world->ymax, current);
                data.func_data.data->push_back({ t, current });
                current = 0;
                ndt = 0;
            }

            m_cSystemPlot.RedrawBuffer();
            m_cSystemPlot.SwapBuffers();
            m_cAnodeCurrentPlot.RedrawBuffer();
            m_cAnodeCurrentPlot.SwapBuffers();
            Invoke([this] () {
                m_cSystemPlot.RedrawWindow();
                m_cAnodeCurrentPlot.RedrawWindow();
            });
        }
    }
    CSimulationDialog::OnSimulation();
}


void CTetrodeDlg::OnBnClickedCheck1()
{
    UpdateData(TRUE);
    data.system_data.point_plot->visible = (m_bPointsVisible == TRUE);
    data.system_data.triangulation_plot->visible = (m_bTriangulationVisible == TRUE);
    data.system_data.dirichlet_cell_plot->visible = (m_bDirichletCellsVisible == TRUE);
    data.isoline_data.plot->visible = (m_bIsolineVisible == TRUE);
}
