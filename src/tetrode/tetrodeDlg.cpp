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
{
    m_hIcon = AfxGetApp()->LoadIcon(IDR_MAINFRAME);
}

void CTetrodeDlg::DoDataExchange(CDataExchange* pDX)
{
    CSimulationDialog::DoDataExchange(pDX);
    DDX_Control(pDX, IDC_MODEL_PLOT, m_cSystemPlot);
    DDX_Control(pDX, IDC_ANODE_CURRENT_OF_TIME_PLOT, m_cITPlot);
    DDX_Control(pDX, IDC_ANODE_CURRENT_OF_POTENTIAL_PLOT, m_cIVPlot);
}

BEGIN_MESSAGE_MAP(CTetrodeDlg, CSimulationDialog)
    ON_WM_PAINT()
    ON_WM_QUERYDRAGICON()
    ON_BN_CLICKED(IDC_BUTTON1, &CTetrodeDlg::OnBnClickedButton1)
    ON_BN_CLICKED(IDC_BUTTON3, &CTetrodeDlg::OnBnClickedButton3)
    ON_BN_CLICKED(IDC_BUTTON2, &CTetrodeDlg::OnBnClickedButton2)
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
    m_eSimulationMode = sm_it;
    UpdateData(TRUE);
    StartSimulationThread();
}


void CTetrodeDlg::OnBnClickedButton2()
{
    StopSimulationThread();
}


void CTetrodeDlg::OnSimulation()
{
    simulation_mode sm;
    Invoke([&] () { sm = m_eSimulationMode; });
    while (m_bWorking)
    {
        Sleep(1000); // stub
    }
    CSimulationDialog::OnSimulation();
}
