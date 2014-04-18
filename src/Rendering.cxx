#include <iostream>

#include <vtkProperty.h>

#include "Rendering.h"

Rendering* Rendering::New() {

  Rendering *r = new Rendering();
 
  return r;

}


void Rendering::update( const std::vector<Vector> &x )
{ 
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> vertices = vtkSmartPointer<vtkCellArray>::New();
  vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
  vtkIdType pid[2];
  for(int i = 0; i < x.size()-1; i++)
  {
    pid[0] = points->InsertNextPoint( x[i][0], x[i][1], x[i][2] );
    pid[1] = points->InsertNextPoint( x[i+1][0], x[i+1][1], x[i+1][2] );
    vertices->InsertNextCell( 1, &pid[0] );
    lines->InsertNextCell( 2, &pid[0] );
  }

  m_poly->SetPoints( points );
  m_poly->SetVerts( vertices );
  m_poly->SetLines( lines );

}

void Rendering::Execute(vtkObject *caller, unsigned long eventId,
                         void * vtkNotUsed(callData))
{
  m_callback();
  m_renderWindow->Render();
}

void Rendering::run()
{
  m_interactor->Start();
}

void Rendering::setup()
{
  // setup vtk pipeline
  m_renderWindow =
    vtkSmartPointer<vtkRenderWindow>::New();

  vtkSmartPointer<vtkRenderer> renderer =
    vtkSmartPointer<vtkRenderer>::New();

  m_renderWindow->AddRenderer(renderer);

  m_interactor =
    vtkSmartPointer<vtkRenderWindowInteractor>::New();

  vtkSmartPointer<vtkInteractorStyle> style = 
    vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();

  m_interactor->SetInteractorStyle(style);

  m_interactor->SetRenderWindow(m_renderWindow);
  m_renderWindow->Render();


  m_poly = vtkSmartPointer<vtkPolyData>::New();
  vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInput(m_poly);


  vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);
  actor->GetProperty()->SetPointSize( 3.0 );

  renderer->AddActor(actor);

  renderer->ResetCamera();
  renderer->GetActiveCamera()->SetPosition(5,0,20);
  renderer->GetActiveCamera()->SetFocalPoint(5,0,0);
  renderer->GetActiveCamera()->Modified();

  // setup timer event
  m_interactor->Initialize();
  m_interactor->AddObserver(vtkCommand::TimerEvent, this);
  m_timerId = m_interactor->CreateRepeatingTimer(20);

}

void Rendering::setCallback( void ( *func)(void) )
{
  m_callback = func;
}



