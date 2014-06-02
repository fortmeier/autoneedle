/**
 * Copyright (C) 2014 Dirk Fortmeier
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE
 *
 */

 #include <iostream>

#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkCellArray.h>
#include <vtkCamera.h>
#include <vtkCommand.h>
#include "mathheader.h"


/**
 * Pretty stupid needle renderer using VTK.
 * A different renderer using for example OGRE could replace it.
 */

class Rendering : public vtkCommand
{
  private:
    virtual void Execute(vtkObject *caller, unsigned long eventId, void * vtkNotUsed(callData));

  public:
    static Rendering *New();
 

    void updateDouble( const std::vector<cml::vector3d> &x );
    void updateFloat( const std::vector<cml::vector3f> &x );

    void run();

    void setup();

    void setCallback( void ( *func)(void) );

  private:
    vtkSmartPointer<vtkRenderWindow> m_renderWindow;
    vtkSmartPointer<vtkPolyData> m_poly;
    vtkSmartPointer<vtkRenderWindowInteractor> m_interactor;
    int m_timerId;
    void (*m_callback)(void);


};


