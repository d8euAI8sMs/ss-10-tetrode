// tetrodeDlg.h : header file
//

#pragma once

#include <util/common/gui/SimulationDialog.h>
#include <util/common/gui/PlotControl.h>

// CTetrodeDlg dialog
class CTetrodeDlg : public CSimulationDialog
{
// Construction
public:
    CTetrodeDlg(CWnd* pParent = NULL);    // standard constructor

// Dialog Data
    enum { IDD = IDD_TETRODE_DIALOG };

    enum simulation_mode { sm_iv, sm_it };

    protected:
    virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support


// Implementation
protected:
    HICON m_hIcon;

    // Generated message map functions
    virtual BOOL OnInitDialog();
    afx_msg void OnPaint();
    afx_msg HCURSOR OnQueryDragIcon();
    DECLARE_MESSAGE_MAP()
    simulation_mode m_eSimulationMode;
public:
    afx_msg void OnBnClickedButton1();
    afx_msg void OnBnClickedButton3();
    afx_msg void OnBnClickedButton2();
    void OnSimulation();
    CPlotControl m_cSystemPlot;
    CPlotControl m_cITPlot;
    CPlotControl m_cIVPlot;
};
