/*
Implementação do plugin Distributed Independent Random Fields SGSIM.

(c) 2014, LPM/UFRGS, Luiz Gustavo Rasera, Péricles Lopes Machado

BASED ON AR2GEMS SGSIM IMPLEMENTATION
*/



#include "report.h"

DPGSGSIM_report::DPGSGSIM_report(QWidget *parent):QDialog(parent){
	ui.setupUi(this);
	self = 0;
    QObject::connect(this,
                     SIGNAL(setText(QString)),
                     this,
                     SLOT(setReport(QString)),
                     Qt::QueuedConnection
                     );

	QObject::connect(this,
		SIGNAL(appendText(QString)),
		this,
		SLOT(appendReport(QString)),
		Qt::QueuedConnection
		);
	this->setAttribute(Qt::WA_DeleteOnClose);

}

DPGSGSIM_report::~DPGSGSIM_report(){

}
