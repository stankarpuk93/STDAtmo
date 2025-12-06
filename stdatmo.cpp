#include <QtCore>
#include <QtGui>
#include <QMessageBox>
#include <cstdlib>
#include <algorithm>

#include "stdatmo.h"
#include "ui_stdatmo.h"
#include "units.h"

// Define a range of constants for
// standard atmosphere calculations
const double Href[8] = {0, 11000, 20000, 32000, 47000, 51000, 71000, 84852}; // geopotential altitude, m
const double Kref[8] = {-0.0065, 0, 0.001, 0.0028, 0, -0.0028, -0.002, -0.002};
const double Tiso[8] = {288.15, 216.66, 216.66, 228.65, 270.65, 270.65, 214.65, 187.65};
const double piso[8] = {101325, 22632, 5474, 868.1, 110.9, 66.94, 3.956, 0.373};
const double g0      = 9.807;
const double R       = 287.04;
const double mu0     = 17.16E-6;
const double S       = 110.6;
const double T0      = 273.15;
const double r0      = 6356.766;
const double gamma   = 1.4;
const double Heps    = 2;       // altitude tolerance useful when units are converted




void compute_speed_of_sound_and_viscosity(double &a, double &mu, double T) {

    // computes speed of sound and dynamic viscosity
    a  = sqrt(gamma * R * T);
    mu = mu0 * pow(T/T0,1.5)*(T0+S)/(T+S);
}

double compute_acceleration_of_gravity(double Hgeo){

    // computes acceleration of gravity wrt altitude
    double g1 = g0 * pow(r0/(r0 + Hgeo / 1000), 2);

    return g1;
}


double compute_geopotential_altitude(double Hgeo){

    // computes a geopotential altitude for a given geometric altitude
    double Hgp = (r0* Hgeo)/(r0 + Hgeo / 1000);

    return Hgp;
}

double compute_geometric_altitude(double Hgp){

    // computes a geometric altitude for a given geopotential altitude
    double Hgeo = (Hgp * r0) / (r0 - Hgp/1000.0);;

    return Hgeo;

}

void append_selected_unit_value(double H, double &g, double &p, double &T, double &rho, double &mu,
                                double &a, QVector<QPair<double, double>> &atmoData,
                                QString pl_type, QString unit_type, QString unit_setup){

    // appends the value of interest into an array to plot
    if (pl_type == "Pressure"){
        if (unit_setup == "SI"){
            if (unit_type == "mmHg"){
                atmoData.append(qMakePair(Units::pascalsToMmHg(p), H));
            }
            else {
                atmoData.append(qMakePair(p, H));
            }
        }
        else {
                atmoData.append(qMakePair(Units::pascalsToMmHg(p), Units::metersToFeet(H)));
            }
        }
    else if (pl_type == "Temperature"){
        if (unit_setup == "SI"){
            if (unit_type == "°C"){
                atmoData.append(qMakePair(Units::kelvinToCelsius(T), H));
            }
            else
            {
                atmoData.append(qMakePair(T, H));
            }
        }
        else{
            if (unit_type == "°F"){
                atmoData.append(qMakePair(Units::kelvinToFahrenheit(T), Units::metersToFeet(H)));
                }
            else{
                atmoData.append(qMakePair(Units::kelvinToRankine(T), Units::metersToFeet(H)));
                }
            }
        }
    else if (pl_type == "Speed of sound"){
        if (unit_setup == "SI"){
            if (unit_type == "km/hr"){
                atmoData.append(qMakePair(Units::metersPerSecondToKilometersPerHour(a), H));
            }
            else{
                atmoData.append(qMakePair(a, H));
            }
        }
        else{
            if (unit_type == "ft/s"){
                atmoData.append(qMakePair(Units::metersPerSecondToFeetPerSecond(a), Units::metersToFeet(H)));
            }
            else if (unit_type == "mi/hr"){
                atmoData.append(qMakePair(Units::metersPerSecondToMilesPerHour(a), Units::metersToFeet(H)));
            }
            else{
                atmoData.append(qMakePair(Units::metersPerSecondToKnots(a), Units::metersToFeet(H)));
            }
            }
        }
        else if (pl_type == "Density"){
            if (unit_setup == "SI"){
                atmoData.append(qMakePair(rho, H));
            }
            else{
                atmoData.append(qMakePair(Units::kgPerM3ToSlugPerFt3(rho), Units::metersToFeet(H)));
            }
        }
    else if (pl_type == "Dynamic viscosity"){
            if (unit_setup == "SI"){
                atmoData.append(qMakePair(mu, H));
            }
            else{
                atmoData.append(qMakePair(Units::pascalSecondToLbfSecondPerFt2(mu), Units::metersToFeet(H)));
            }
        }
    else{
            if (unit_setup == "SI"){
                atmoData.append(qMakePair(g, H));
            }
            else{
                atmoData.append(qMakePair(Units::metersPerSecondToFeetPerSecond(g), Units::metersToFeet(H)));
            }
        }


}


void compute_standard_atmosphere(double H, double dISA, double &g, double &p, double &T, double &rho, double &mu, double &a)
{

    // define all necessary constants
    double rho_iso;
    int Htest_ind = 1;

    // find the index of a required altitude
    while (H > Href[Htest_ind])
    {
        Htest_ind += 1;
    }
    Htest_ind -= 1;

    double K  = Kref[Htest_ind];
    double H1 = Href[Htest_ind];
    double T1 = Tiso[Htest_ind];

    if (K == 0)
    {
        // isothermal region
        double alpha = exp(-g0*(H - H1) / R / T1);
        T   = T1 + dISA;
        p   = piso[Htest_ind] * alpha;
        rho_iso = piso[Htest_ind] / R / T;
        rho = rho_iso * alpha;

    }
    else
    {
        // gradient region
        double alpha = -g0 / R / K;
        double T_noISA = T1 + K * (H - H1);
        T = T_noISA + dISA;
        p = piso[Htest_ind] * pow(T_noISA / T1, alpha);
        rho = p / R / T;
    }

    compute_speed_of_sound_and_viscosity(a, mu, T);
    g = compute_acceleration_of_gravity(H);
}


STDAtmo::STDAtmo(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::STDAtmo)
{
    ui->setupUi(this);


    // plot default parameters for the plotting widget
    ui->PlottingWidget->xAxis->setLabel("Temperature, K");
    ui->PlottingWidget->yAxis->setLabel("Altitude, m");
    ui->PlottingWidget->xAxis->setRange(0,288);
    ui->PlottingWidget->yAxis->setRange(0,86000);
    ui->PlottingWidget->replot();


    // write an information section
    ui -> infoBrowser -> setOpenExternalLinks(true);
    QString content = R"(
    <h1>STDAtmo</h1>

    <h2>Description</h2>
    <p>The STDAtmo computes standard atmosphere properties based on the the 1976 model.
       All typical air properties with ISA deviations can be computed up to an altitude of 86 kilometers.</p>

    <h2>Capabilities</h2>
    <p>STDAtmo has two features:
    <ul>
    <li> Calculations of standard atmosphere properties at a given altitude and temperature deviations in SI and US units.
    <li> Plotting of standard atmosphere properties in a user-defined range and prescribed units.
    </ul>
    </p>

    <h2>References</h2>
    <ul>
    <li><a href="https://ntrs.nasa.gov/api/citations/19770009539/downloads/19770009539.pdf">U.S. Standard Atmosphere 1976 (NASA-TM-X-74335)</a> - standard atmosphere methodology</li>
    </ul>

    <h3>Author</h3>
    <p>  <li><a href="https://www.linkedin.com/in/stankarpuk/">Stanislav Karpuk</a></li></p>
    )";

    ui -> infoBrowser -> setHtml(content);

    // initialize the initial plot and the default range
    set_default_plot_inputs();

}

STDAtmo::~STDAtmo()
{
    delete ui;

}


void STDAtmo::on_Reset_clicked()
{
    computed = false;

    // resets all units
    QList<QComboBox*> comboBoxes = ui->centralwidget->findChildren<QComboBox*>();
    for (QComboBox* comboBox : std::as_const(comboBoxes)) {
        comboBox->setCurrentIndex(0);  // Set to first item
    }

    // clears lineEdits
    QList<QLineEdit*> lineEdits = ui->centralwidget->findChildren<QLineEdit*>();
    for (QLineEdit* lineEdit : lineEdits) {
        lineEdit->clear();
    }

    // set plotting inputs to default
    set_default_plot_inputs();

    // clean the plotting widget

    ui->PlottingWidget->xAxis->setLabel("Temperature");
    ui->PlottingWidget->yAxis->setLabel("Altitude");
    ui->PlottingWidget->replot();

}


void STDAtmo::set_default_plot_inputs(){

    // resets the default set of inputs for plotting
    // (100 points by default from 0 to 85 km)
    ui -> AltitudeMinInp -> setText(QString::number(0));
    ui -> AltitudeMaxInp -> setText(QString::number(86000));
    ui -> DAltInp -> setText(QString::number(1000));
}

void STDAtmo::on_Compute_clicked()
{

    // Read the input values and convert to double
    QString dISAtext = ui -> dISAInp -> text();
    QString Alttext = ui -> AltitudeInp -> text();

    dISA = dISAtext.toDouble();
    H  = Alttext.toDouble();            // geometric altitude

    // compute a geopotential altitude
    Hgp = compute_geopotential_altitude(H);

    // Read the units and convert to SI
    QString AltUnits = ui -> AltUnits -> currentText();
    QString dISAUnits = ui -> TempUnits -> currentText();

    if (AltUnits == "ft") { Hgp /= 3.28084; }
    if (dISAUnits == "°F / °R") { dISA = dISA * 5/9; }

    // Check the altitude limit
    if (Hgp > Href[7]){
        if (Hgp > Href[7]+ Heps)     // adds a 0.5 m tolerance
        {
            QMessageBox::information(this,"Title","You are aiming too high, son...\n Enter an altitude below 86 km");
        }
        else{
            Hgp = Href[7];
        }
    }


    compute_standard_atmosphere(Hgp, dISA, g, p, T, ro, mu, a);
    computed = true;

    // output results in appropriate units
    output_results();

}

void STDAtmo::on_GopAltCombo_currentIndexChanged(){output_results();}

void STDAtmo::on_PressureCombo_currentIndexChanged(){output_results();}

void STDAtmo::on_TemperatureCombo_currentIndexChanged(){output_results();}

void STDAtmo::on_DensityCombo_currentIndexChanged(){output_results();}

void STDAtmo::on_SpeedSoundCombo_currentIndexChanged(){output_results();}

void STDAtmo::on_ViscosityCombo_currentIndexChanged(){output_results();}

void STDAtmo::on_GravityCombo_currentIndexChanged(){output_results();}

void STDAtmo::output_results(){

    /*
     The code output the data based on the combobox and
        local units indices
    */

    if (!computed) return;

    int selectedText;

    // Altitude output
    selectedText = ui-> GopAltCombo -> currentIndex();
    switch (selectedText){
    case 0: ui -> AltitudeOutput -> setText(QString::number(Hgp)); break;
    case 1: ui -> AltitudeOutput -> setText(QString::number(Units::metersToFeet(Hgp))); break;
    }

    // Pressure output
    selectedText = ui->PressureCombo->currentIndex();
    switch (selectedText) {
        case 0: ui -> PressureOutput -> setText(QString::number(p));    break;
        case 1: ui -> PressureOutput -> setText(QString::number(Units::pascalsToPsi(p)));   break;
        case 2: ui -> PressureOutput -> setText(QString::number(Units::pascalsToMmHg(p)));   break;
    }

    // Temperature output
    selectedText = ui->TemperatureCombo->currentIndex();
    switch (selectedText) {
        case 0: ui -> TemperatureOutput -> setText(QString::number(T)); break;
        case 1: ui -> TemperatureOutput -> setText(QString::number(Units::kelvinToCelsius(T))); break;
        case 2: ui -> TemperatureOutput -> setText(QString::number(Units::kelvinToFahrenheit(T)));  break;
        case 3: ui -> TemperatureOutput -> setText(QString::number(Units::kelvinToRankine(T))); break;
    }

    // Density output
    selectedText = ui->DensityCombo->currentIndex();
    switch (selectedText) {
        case 0: ui -> DensityOutput -> setText(QString::number(ro));    break;
        case 1: ui -> DensityOutput -> setText(QString::number(Units::kgPerM3ToSlugPerFt3(ro)));    break;
    }

    // Speed of sound
    selectedText = ui->SpeedSoundCombo->currentIndex();
    switch (selectedText) {
        case 0: ui -> SoundSpeedOutput -> setText(QString::number(a));  break;
        case 1: ui -> SoundSpeedOutput -> setText(QString::number(Units::metersPerSecondToFeetPerSecond(a)));   break;
        case 2: ui -> SoundSpeedOutput -> setText(QString::number(Units::metersPerSecondToKilometersPerHour(a)));   break;
        case 3: ui -> SoundSpeedOutput -> setText(QString::number(Units::metersPerSecondToMilesPerHour(a)));    break;
        case 4: ui -> SoundSpeedOutput -> setText(QString::number(Units::metersPerSecondToKnots(a)));   break;
    }

    // Dynamic viscosoty
    selectedText = ui->ViscosityCombo->currentIndex();
    switch (selectedText) {
        case 0: ui -> ViscosityOutput -> setText(QString::number(mu));  break;
        case 1: ui -> ViscosityOutput -> setText(QString::number(Units::pascalSecondToLbfSecondPerFt2(mu)));    break;
    }

    // Acceleration of gravity
    selectedText = ui->GravityCombo->currentIndex();
    switch (selectedText) {
    case 0: ui -> GravityOutput -> setText(QString::number(g));  break;
    case 1: ui -> GravityOutput -> setText(QString::number(Units::metersPerSecondToFeetPerSecond(g)));    break;
    }
}

void STDAtmo::on_actionExit_triggered()
{
    QApplication::quit();
}


void STDAtmo::on_UnitTypePlot_currentIndexChanged(){

    // defines units for plotting
    update_plotting_units();
}


void STDAtmo::on_UnitSetup_currentIndexChanged(){

    // defines units for plotting
    update_plotting_units();

    // updates input units and EditText inputs according to specofied units
    QString unit;

    double Harr[3];
    QString HminText = ui -> AltitudeMinInp -> text();
    QString HmaxText = ui -> AltitudeMaxInp -> text();
    QString dHText = ui -> DAltInp -> text();
    Harr[0] = HminText.toDouble();
    Harr[1] = HmaxText.toDouble();
    Harr[2] = dHText.toDouble();

    int UnitsIndex = ui -> UnitSetup -> currentIndex();
    if (UnitsIndex == 1){
        unit = "ft";
        for (int i = 0; i < 3; ++i) {
            Harr[i] = Units::metersToFeet(Harr[i]);
        }
    }
    else {
        unit = "m";
        for (int i = 0; i < 3; ++i) {
            Harr[i] = Units::feetToMeter(Harr[i]);
        }
    }
    // write converted values
    ui -> AltitudeMinInp -> setText(QString::number(Harr[0]));
    ui -> AltitudeMaxInp -> setText(QString::number(Harr[1]));
    ui -> DAltInp -> setText(QString::number(Harr[2]));

    // Assign new labels
    ui -> minAltUnit -> setText(unit);
    ui -> maxAltUnit -> setText(unit);
    ui -> dAltUnit -> setText(unit);

}



void STDAtmo::update_plotting_units(){

    QString unitText = ui -> UnitSetup -> currentText();
    QString comboText = ui -> UnitTypePlot -> currentText();

    //Modify selected units
    ui -> UnitPlot -> clear();
    if (comboText == "Pressure"){
        if (unitText == "SI"){
            QStringList list=(QStringList()<<"Pa"<<"mmHg");
            ui -> UnitPlot -> addItems(list);
        }
        else {
            QStringList list=(QStringList()<<"psi");
            ui -> UnitPlot -> addItems(list);
        }
    }
    else if (comboText == "Temperature"){
        if (unitText == "SI"){
            QStringList list=(QStringList()<<"K"<<"°C");
            ui -> UnitPlot -> addItems(list);
        }
        else {
            QStringList list=(QStringList()<<"°F"<<"°R");
            ui -> UnitPlot -> addItems(list);
        }
    }
    else if (comboText == "Density"){
        if (unitText == "SI"){
            QStringList list=(QStringList()<<"kg/m3");
            ui -> UnitPlot -> addItems(list);
        }
        else {
            QStringList list=(QStringList()<<"lb/ft3");
            ui -> UnitPlot -> addItems(list);
        }
    }
    else if (comboText == "Speed of sound"){
        if (unitText == "SI"){
            QStringList list=(QStringList()<<"m/s"<<"km/hr");
            ui -> UnitPlot -> addItems(list);
        }
        else {
            QStringList list=(QStringList()<<"ft/s"<<"mi/hr"<<"knots");
            ui -> UnitPlot -> addItems(list);
        }
    }
    else if (comboText == "Dynamic viscosity"){
        if (unitText == "SI"){
            QStringList list=(QStringList()<<"Pa.s");
            ui -> UnitPlot -> addItems(list);
        }
        else {
            QStringList list=(QStringList()<<"lbf.s/ft2");
            ui -> UnitPlot -> addItems(list);
        }
    }
    else if (comboText == "Acceleration of gravity"){
        if (unitText == "SI"){
            QStringList list=(QStringList()<<"m/s2");
            ui -> UnitPlot -> addItems(list);
        }
        else {
            QStringList list=(QStringList()<<"ft/s2");
            ui -> UnitPlot -> addItems(list);
        }
    }

}


void STDAtmo::on_Plot_graph_clicked()
{

    ui->PlottingWidget->clearGraphs();

    // plots standard atmosphere properties
    QString HminText = ui -> AltitudeMinInp -> text();
    QString HmaxText = ui -> AltitudeMaxInp -> text();
    QString dHText = ui -> DAltInp -> text();
    QString pl_type = ui -> UnitTypePlot -> currentText();
    QString unit_inp = ui -> UnitPlot -> currentText();
    QString unit_def = ui -> UnitSetup -> currentText();

    hbegin = HminText.toDouble();
    hend   = HmaxText.toDouble();
    dh     = dHText.toDouble();

    // Configure plot
    QString yLabel = (unit_def == "SI") ? "Altitude, m" : "Altitude, ft";
    ui->PlottingWidget->yAxis->setLabel(yLabel);
    ui->PlottingWidget->xAxis->setLabel(pl_type + ", " + unit_inp);

    // Convert input units to SI for calculations
    if (unit_def == "US") {
        // Convert from feet to meters for internal calculations
        hbegin = Units::feetToMeter(hbegin);
        hend = Units::feetToMeter(hend);
        dh = Units::feetToMeter(dh);
    }

    // compute geopotential altitudes
    double hbeg_geo = compute_geopotential_altitude(hbegin);
    double hend_geo = compute_geopotential_altitude(hend);

    QVector<QPair<double, double>> atmoData;

    int Npoints = static_cast<int>(std::ceil((hend - hbegin) / dh)) + 1;        // number of points to plot except for the lazers boundaries
    atmoData.reserve(Npoints + 8);   // reserves extra 8 points for potential layers

    // compute points to plot
    for (double X = hbegin; X <= hend + 1e-6; X += dh) {
        double Hgp = compute_geopotential_altitude(X);
        compute_standard_atmosphere(Hgp, dISA, g, p, T, ro, mu, a);
        append_selected_unit_value(compute_geometric_altitude(Hgp), g, p, T, ro, mu, a, atmoData,
                                   pl_type, unit_inp, unit_def);
    }

    // determine how many layers and points are required
    for (int i = 0; i < 8; ++i) {
        if (Href[i] > hbeg_geo && Href[i] < hend_geo) {
            compute_standard_atmosphere(Href[i], dISA, g, p, T, ro, mu, a);
            append_selected_unit_value(compute_geometric_altitude(Href[i]), g, p, T, ro, mu, a, atmoData,
                                       pl_type, unit_inp, unit_def);
        }
    }

    // sort the data
    std::sort(atmoData.begin(), atmoData.end(),
              [](const QPair<double, double>& a, const QPair<double, double>& b) {
                  return a.second < b.second;
              });


    // plot results segment-by-segment
    QVector<double> Xmax, Ymax, Xmin, Ymin;
    int datacount = 0, gp_alt;
    int data_size = atmoData.size();
    for (int i = 1; i < 8; i++){
        QVector<double> xData, yData;
        xData.reserve(data_size);
        yData.reserve(data_size);
        if (datacount < data_size) {
            do {
                // compute geopotential altitude depending on units specified
                if (unit_def == "US") {
                    gp_alt = compute_geopotential_altitude(Units::feetToMeter(atmoData[datacount].second));
                }
                else{
                    gp_alt = compute_geopotential_altitude(atmoData[datacount].second);
                }
                // append calculated values into separate vectors
                xData.append(atmoData[datacount].second);
                yData.append(atmoData[datacount].first);
                datacount += 1;
            } while (gp_alt + Heps < Href[i] && datacount < data_size);         // used tolerance to avoid issue while converting units

            // plot the graph
            if (!xData.isEmpty()) {
                ui->PlottingWidget->addGraph();
                ui->PlottingWidget->graph(i-1)->setData(yData, xData);

                Xmax.append(*std::max_element(yData.begin(), yData.end()));
                Ymax.append(*std::max_element(xData.begin(), xData.end()));
                Xmin.append(*std::min_element(yData.begin(), yData.end()));
                Ymin.append(*std::min_element(xData.begin(), xData.end()));

                ui->PlottingWidget-> xAxis -> setRange(*std::min_element(Xmin.begin(), Xmin.end()) * 0.99,
                                            *std::max_element(Xmax.begin(), Xmax.end()) * 1.01);
                ui->PlottingWidget-> yAxis -> setRange(*std::min_element(Ymin.begin(), Ymin.end()),
                                            *std::max_element(Ymax.begin(), Ymax.end()) * 1.01);

                ui->PlottingWidget->replot();
            }
            datacount -= 1;
        }
        else break;

    }

}


