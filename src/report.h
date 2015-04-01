/*
Implementação do plugin Distributed Independent Random Fields SGSIM.

(c) 2014, LPM/UFRGS, Luiz Gustavo Rasera, Péricles Lopes Machado

BASED ON AR2GEMS SGSIM IMPLEMENTATION
*/


#ifndef DPGSGSIM_REPORT_H_
#define DPGSGSIM_REPORT_H_

#include "common.h"
#include <QDialog>
#include <QMainWindow>
#include <QCloseEvent>
#include <ui_report.h>


class DPGSGSIM_report : public QDialog
{
	Q_OBJECT

public:
    DPGSGSIM_report(QWidget *parent = 0);
    ~DPGSGSIM_report();

    void write(QString msg){
        emit setText(msg);
    }

	void append(QString msg) {
		emit appendText(msg);
	}

signals:
    void setText(QString text);
	void appendText(QString text);

public slots:

    void setReport(QString text) {
        this->ui.report_text->setText(text);
        this->show();
    }

	void appendReport(QString text) {
		this->ui.report_text->append(text);
		this->show();
	}

    void setSelf(DPGSGSIM_report* r) {
		self = r;
	}

	void autoDestroy() {
		if (self) delete self;
	}

	void closeEvent(QCloseEvent *event) {
		if (self) delete self;
	}

private:
    DPGSGSIM_report* self;
    Ui::DPGSGSIM_report ui;
};

#endif /* LPM_DPGSGSIM_REPORT_H_ */
