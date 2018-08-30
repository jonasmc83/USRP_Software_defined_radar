/*
    USRP_Software_defined_radar is a software for real time sampling, processing, display and storing
    Copyright (C) 2018  Jonas Myhre Christiansen <jonas-myhre.christiansen@ffi.no>
	
    This file is part of USRP_Software_defined_radar.

    USRP_Software_defined_radar is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    USRP_Software_defined_radar is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with USRP_Software_defined_radar.  If not, see <https://www.gnu.org/licenses/>.
*/

#include "mainwindow.h"
#include <QApplication>
#include <QDebug>
#include <QMessageBox>
#include <exception>

class MyApplication : public QApplication
{
//    Q_OBJECT
public:
    MyApplication(int &argc, char **argv) : QApplication(argc, argv) {}

    bool notify(QObject* receiver, QEvent* event) //Q_DECL_OVERRIDE
    {
        try {
            return QApplication::notify(receiver, event);
        } catch (std::exception &e) {
            // Handle the desired exception type
            QMessageBox::critical(this->activeWindow(), "Exceptions thrown...", e.what());
            this->exit(-1);
        } catch (...) {
            // Handle the rest
        }

         return false;
     }
};

int UHD_SAFE_MAIN(int argc, char *argv[])
{
	std::cout << "USRP_Software_defined_radar  Copyright (C) 2018  Jonas Myhre Christiansen <jonas-myhre.christiansen@ffi.no>" << std::endl;
    std::cout << "This program comes with ABSOLUTELY NO WARRANTY." << std::endl;
    std::cout << "This is free software, and you are welcome to redistribute it" << std::endl;
    std::cout << "under certain conditions." << std::endl;
	
    MyApplication a(argc, argv);
    MainWindow w;
    w.show();

    return a.exec();
}
