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



class Rendering : public vtkCommand
{
  public:
    static Rendering *New();
 
    virtual void Execute(vtkObject *caller, unsigned long eventId,
                         void * vtkNotUsed(callData));

    void update( const std::vector<Vector> &x );

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


