#include "stdatmo.h"

#include <QIcon>
#include <QApplication>
#include <QDebug>    // Для отладки
#include <QFile>     // Для проверки существования файла
#include <QDir>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);

    // set the icon up
    a.setWindowIcon(QIcon(":/STDAtmo_logo.png"));

    // get the application data
    QApplication::setApplicationName("STDAtmo");
    QApplication::setApplicationVersion("1.0.0");
    QApplication::setOrganizationName("Stanislav Karpuk");
    QApplication::setOrganizationDomain("stankarpuk93@gmail.com");


    STDAtmo w;
    w.show();
    return a.exec();
}
