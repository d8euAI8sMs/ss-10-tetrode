// tetrode.h : main header file for the PROJECT_NAME application
//

#pragma once

#ifndef __AFXWIN_H__
    #error "include 'stdafx.h' before including this file for PCH"
#endif

#include "resource.h"        // main symbols


// CTetrodeApp:
// See tetrode.cpp for the implementation of this class
//

class CTetrodeApp : public CWinApp
{
public:
    CTetrodeApp();

// Overrides
public:
    virtual BOOL InitInstance();

// Implementation

    DECLARE_MESSAGE_MAP()
};

extern CTetrodeApp theApp;
