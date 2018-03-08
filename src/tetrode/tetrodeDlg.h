// tetrodeDlg.h : header file
//

#pragma once

#include <util/common/gui/SimulationDialog.h>

// CTetrodeDlg dialog
class CTetrodeDlg : public CSimulationDialog
{
// Construction
public:
    CTetrodeDlg(CWnd* pParent = NULL);    // standard constructor

// Dialog Data
    enum { IDD = IDD_TETRODE_DIALOG };

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
};
