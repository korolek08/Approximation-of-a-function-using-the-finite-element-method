#include <QtWidgets/QApplication>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QMessageBox>
#include <QtWidgets/QAction>
#include "window.h"

int main (int argc, char *argv[]) {
    
    QApplication app(argc, argv);
    QLocale::setDefault(QLocale::C);
    QMainWindow *window = new QMainWindow;
    Window *graph_area = new Window(window);
    QMenuBar *tool_bar = new QMenuBar(window);
    QAction *action;
     
    if (graph_area->commandLine(argc, argv)) {
        printf("Usage: %s a b c d nx ny mx my k eps maxit p", argv[0]);
        //app.exec ();
        //graph_area->delete_mem();
        //delete graph_area;
        delete tool_bar;
        //delete window;
        return -1;          
    }
    
    int p = graph_area->return_p();
    Args* ar = new Args[p];
    for (int k = 0; k < p; k++) {
        ar[k].arguments = argv[0];
        ar[k].n_thread = k;
        ar[k].p = p;
    }
    //printf("%d n\n", graph_area->n);
    graph_area->ar = ar;
    graph_area->init_r();
    graph_area->init_mem();
    graph_area->connection_continuous();
    graph_area->add_param(ar);
    //printf("r\n");
    //printf("rrrr\n");
    action = tool_bar->addAction ("&Change function", graph_area, SLOT(change_func()));
    action->setShortcut(QString("0"));
    action = tool_bar->addAction("&Change state", graph_area, SLOT(change_state()));
    action->setShortcut(QString("1"));
    action = tool_bar->addAction("&Up scale", graph_area, SLOT(up_scale()));
    action->setShortcut(QString("2"));
    action = tool_bar->addAction("&Down scale", graph_area, SLOT(down_scale()));
    action->setShortcut(QString("3"));
    action = tool_bar->addAction("&Up n", graph_area, SLOT(up_n()));
    action->setShortcut(QString("4"));
    action = tool_bar->addAction("&Down n", graph_area, SLOT(down_n()));
    action->setShortcut(QString("5"));
    action = tool_bar->addAction("&Up perturbation", graph_area, SLOT(up_perturbation()));
    action->setShortcut(QString("6"));
    action = tool_bar->addAction("&Down perturbation", graph_area, SLOT(down_perturbation()));
    action->setShortcut(QString("7"));
    action = tool_bar->addAction("&Up visualization", graph_area, SLOT(up_visualization()));
    action->setShortcut(QString("8"));
    action = tool_bar->addAction("&Down visualization", graph_area, SLOT(down_visualization()));
    action->setShortcut(QString("9"));
    action = tool_bar->addAction("E&xit", window, SLOT(close()));
    action->setShortcut(QString("Ctrl+X"));

    tool_bar->setMaximumHeight (30);

    window->setMenuBar (tool_bar);
    window->setCentralWidget (graph_area);
    window->setWindowTitle ("Graph");

    for (int kk=0; kk<p; kk++) {
        if (pthread_create(&ar[kk].tid, 0, continue_thread, ar+kk)) {
            printf("Cannot create thread %d\n", kk);
            abort();
        }
        //printf("TO RESTART\n");
    }
    
    graph_area->restart();
    //ar[0].tid=pthread_self();
    //thread_func(ar+0);
    
    window->show();


    app.exec();
    delete action;
    //graph_area->delete_mem();
    delete graph_area;
    delete tool_bar;
    delete window;
    delete [] ar;
    return 0;
}
