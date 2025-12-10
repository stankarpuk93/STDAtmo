#ifndef STDATMO_H
#define STDATMO_H

#include <QMainWindow>
#include <QVector>
#include <QMainWindow>
#include <QComboBox>
#include <QLineEdit>

QT_BEGIN_NAMESPACE
namespace Ui {
class STDAtmo;
}
QT_END_NAMESPACE

class STDAtmo : public QMainWindow
{
    Q_OBJECT

public:
    STDAtmo(QWidget *parent = nullptr);
    ~STDAtmo();

private slots:
    void on_Reset_clicked();

    void on_Compute_clicked();

    void on_PressureCombo_currentIndexChanged();

    void on_TemperatureCombo_currentIndexChanged();

    void on_DensityCombo_currentIndexChanged();

    void on_SpeedSoundCombo_currentIndexChanged();

    void on_ViscosityCombo_currentIndexChanged();

    void output_STDAresults();

    void on_actionExit_triggered();

    void on_UnitTypePlot_currentIndexChanged();

    void update_plotting_units();

    void on_UnitSetup_currentIndexChanged();

    void set_default_plot_inputs();

    void on_GopAltCombo_currentIndexChanged();

    void on_Plot_graph_clicked();

    void on_GravityCombo_currentIndexChanged();

    void on_Reset_Yp_clicked();

    void on_AltUnits_Yp_2_currentIndexChanged();

    void on_Compute_Yp_clicked();

    void convert_input_values_Ypp(double &Uinf, double &rhoinf, double &muinf, double &Lref);

    void on_Reset_SP_clicked();

    void on_Compute_SP_clicked();

    void convert_airspeed(QComboBox* comboBox, QLineEdit* outputEdit, double airspeed);

    double convert_airspeed_Input(QComboBox* comboBox, double airspeed);

    void on_CASradioButton_clicked();

    void on_MachradioButton_clicked();

    void on_EASradioButton_clicked();

    void on_TASradioButton_clicked();

    void on_TASCombo_currentIndexChanged();

    void on_EASCombo_currentIndexChanged();

    void on_CASCombo_currentIndexChanged();

private:
    Ui::STDAtmo *ui;

    // Initialize values to be computed
    double p, T, ro, mu, a, Hgp, g, ds, VTAS, VCAS, VEAS, Mach;
    double H = 0, dISA = 0;
    bool computed = false;
    double hbegin, hend, dh;
    int Npoints;

};
#endif // STDATMO_H
