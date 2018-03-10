// tetrodeDlg.h : header file
//

#pragma once

#include <util/common/gui/SimulationDialog.h>
#include <util/common/gui/PlotControl.h>

#include "model.h"

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
    model::model_data data;
public:
    afx_msg void OnBnClickedButton1();
    afx_msg void OnBnClickedButton3();
    afx_msg void OnBnClickedButton2();
    void OnSimulation();
    CPlotControl m_cSystemPlot;
    CPlotControl m_cAnodeCurrentPlot;
    afx_msg void OnBnClickedCheck1();
    BOOL m_bPointsVisible;
    BOOL m_bTriangulationVisible;
    BOOL m_bDirichletCellsVisible;
    BOOL m_bIsolineVisible;
    UINT m_nIsolineCount;
    UINT m_nAnodePotentialSamples;
    double m_lfIsolineDelta;
    double m_lfAccuracyGoal;
    double m_lfAnodeBeginPotential;
    double m_lfAnodeEndPotential;
};
