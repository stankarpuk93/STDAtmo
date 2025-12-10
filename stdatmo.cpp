#include <QtCore>
#include <QtGui>
#include <QMessageBox>
#include <cstdlib>
#include <algorithm>
#include <array>
#include <fstream>

#include "stdatmo.h"
#include "ui_stdatmo.h"
#include "units.h"
#include "standardAtmosphere.h"


// Define necessary constants and sole global variables
const double Heps0  = 2.0;
QVector<QPair<double, double>> atmoData;
QString plot_unit_type, plot_unitsX, plot_unitsY;



bool isDouble(const QString& str) {

    // the function checks if the inputs is a number

    bool ok;
    str.toDouble(&ok);
    return ok;
}


void append_selected_unit_value(double H, double &g, double &p, double &T, double &rho, double &mu,
                                double &a, QVector<QPair<double, double>> &atmoData,
                                QString pl_type, QString unit_type, QString unit_setup){

    // appends the value of interest into an array to plot
    if (pl_type == "Pressure"){
        plot_unit_type = "Pressure";
        if (unit_setup == "SI"){
            plot_unitsY = "m";
            if (unit_type == "mmHg"){
                plot_unitsX = "mmHg";
                atmoData.append(qMakePair(Units::pascalsToMmHg(p), H));
            }
            else {
                plot_unitsX = "Pa";
                atmoData.append(qMakePair(p, H));
            }
        }
        else {
            plot_unitsY = "ft";
            plot_unitsX = "lb/ft2";
            atmoData.append(qMakePair(Units::pascalSecondToLbfSecondPerFt2(p), Units::metersToFeet(H)));
            }
        }
    else if (pl_type == "Temperature"){
        plot_unit_type = "Temperature";
        if (unit_setup == "SI"){
            plot_unitsY = "m";
            if (unit_type == "°C"){
                plot_unitsX = "°C";
                atmoData.append(qMakePair(Units::kelvinToCelsius(T), H));
            }
            else
            {
                plot_unitsX = "K";
                atmoData.append(qMakePair(T, H));
            }
        }
        else{
            plot_unitsY = "m";
            if (unit_type == "°F"){
                plot_unitsX = "°F";
                atmoData.append(qMakePair(Units::kelvinToFahrenheit(T), Units::metersToFeet(H)));
                }
            else{
                plot_unitsX = "°R";
                atmoData.append(qMakePair(Units::kelvinToRankine(T), Units::metersToFeet(H)));
                }
            }
        }
    else if (pl_type == "Speed of sound"){
        plot_unit_type = "Speed of sound";
        if (unit_setup == "SI"){
            plot_unitsY = "m";
            if (unit_type == "km/hr"){
                plot_unitsX = "km/hr";
                atmoData.append(qMakePair(Units::metersPerSecondToKilometersPerHour(a), H));
            }
            else{
                plot_unitsX = "m/s";
                atmoData.append(qMakePair(a, H));
            }
        }
        else{
            plot_unitsY = "ft";
            if (unit_type == "ft/s"){
                plot_unitsX = "ft/s";
                atmoData.append(qMakePair(Units::metersPerSecondToFeetPerSecond(a), Units::metersToFeet(H)));
            }
            else if (unit_type == "mi/hr"){
                plot_unitsX = "mi/hr";
                atmoData.append(qMakePair(Units::metersPerSecondToMilesPerHour(a), Units::metersToFeet(H)));
            }
            else{
                plot_unitsX = "knots";
                atmoData.append(qMakePair(Units::metersPerSecondToKnots(a), Units::metersToFeet(H)));
            }
            }
        }
        else if (pl_type == "Density"){
            plot_unit_type = "Density";
            if (unit_setup == "SI"){
                plot_unitsY = "m";
                plot_unitsX = "kg/m3";
                atmoData.append(qMakePair(rho, H));
            }
            else{
                plot_unitsY = "ft";
                plot_unitsX = "slug/ft3";
                atmoData.append(qMakePair(Units::kgPerM3ToSlugPerFt3(rho), Units::metersToFeet(H)));
            }
        }
    else if (pl_type == "Dynamic viscosity"){
            plot_unit_type = "Dynamic viscosity";
            if (unit_setup == "SI"){
                plot_unitsY = "m";
                plot_unitsX = "Pa.s";
                atmoData.append(qMakePair(mu, H));
            }
            else{
                plot_unitsY = "ft";
                plot_unitsX = "lb.s/ft2";
                atmoData.append(qMakePair(Units::pascalSecondToLbfSecondPerFt2(mu), Units::metersToFeet(H)));
            }
        }
    else{
            if (unit_setup == "SI"){
                plot_unitsY = "m";
                plot_unitsX = "m/s2";
                atmoData.append(qMakePair(g, H));
            }
            else{
                plot_unitsY = "ft";
                plot_unitsX = "ft/s2";
                atmoData.append(qMakePair(Units::metersPerSecondToFeetPerSecond(g), Units::metersToFeet(H)));
            }
        }
}


STDAtmo::STDAtmo(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::STDAtmo)
{
    ui->setupUi(this);


    // plot default parameters for the plotting widget
    ui->PlottingWidget->xAxis->setLabel(" ");
    ui->PlottingWidget->yAxis->setLabel("Altitude, m");
    ui->PlottingWidget->xAxis->setRange(0,1);
    ui->PlottingWidget->yAxis->setRange(0,86000);
    ui->PlottingWidget->replot();


    // write infor for the Y+ calculator
    QString footnote = R"(
    1. Y+ calculator is used only for pinitial first mesh layer estimations <br>
    2. Refer to the documentation for the methodology description
    )";

    ui -> infoBrowser_Yp -> setHtml(footnote);


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
    <p><li><a href="https://www.linkedin.com/in/stankarpuk/">Stanislav Karpuk</a></li></p>
    )";

    ui -> infoBrowser -> setHtml(content);

    // initialize the calculation description for the airspeeds tab
    QString instructions = R"(
    <strong>Instructions</strong><br>
    1. Input the altitude and temperature deviations <br>
    2. Select the airspeed to convert from by selecting a radio button <br>
    3. Enter the airspeed of interest and click "Compute" <br><br>
    * Refer to the documentation for the methodology description
    )";
    ui -> infoBrowser_SP -> setHtml(instructions);

    QString STDAref = R"(
    1. Refer to the documentation for the methodology description
    )";
    ui -> infoBrowser_STD -> setHtml(STDAref);


    // initialize the initial plot and the default range
    set_default_plot_inputs();

    // set default airspeed radio button to CAS
    ui -> CASradioButton ->click();

    // disable all velocity inputs except for the CAS
    ui->TASInpOut->setReadOnly(true);
    ui->TASInpOut->setStyleSheet("background-color: #222222;");
    ui->EASInpOut->setReadOnly(true);
    ui->EASInpOut->setStyleSheet("background-color: #222222;");
    ui->MachInpOut->setReadOnly(true);
    ui->MachInpOut->setStyleSheet("background-color: #222222;");
}

STDAtmo::~STDAtmo()
{
    delete ui;

}


void STDAtmo::checkAirspeed(double airspeed){

    if (airspeed < 0 ) {
        QMessageBox::warning(this, "Error", "Enter a positive airspeed");
        return;
    }
}

double STDAtmo::convert_airspeed_Input(QComboBox* comboBox, double airspeed) {

    // reads the combobox index and converts the airspeed to m/s

    int selectedIndex = comboBox->currentIndex();  // Use parameter
    double result;

    switch (selectedIndex) {
    case 0: result = airspeed; break;
    case 1: result = Units::feetPerSecondToMetersPerSecond(airspeed); break;
    case 2: result = Units::kilometersPerHourToMetersPerSecond(airspeed); break;
    case 3: result = Units::milesPerHourToMetersPerSecond(airspeed); break;
    case 4: result = Units::KnotsToMetersPerSecond(airspeed); break;
    }

    return result;
}

void STDAtmo::convert_airspeed(QComboBox* comboBox, QLineEdit* outputEdit, double airspeed) {

    int selectedIndex = comboBox->currentIndex();  // Use parameter
    double result;

    QString init_text = outputEdit -> text();
    if (init_text.isEmpty() && !computed) return;
    else{
        switch (selectedIndex) {
        case 0: result = airspeed; break;
        case 1: result = Units::metersPerSecondToFeetPerSecond(airspeed); break;
        case 2: result = Units::metersPerSecondToKilometersPerHour(airspeed); break;
        case 3: result = Units::metersPerSecondToMilesPerHour(airspeed); break;
        case 4: result = Units::metersPerSecondToKnots(airspeed); break;
        }

        outputEdit->setText(QString::number(result, 'g', 4));
    }

}


void STDAtmo::on_Reset_clicked()
{

    // resets all units and clears all inputs nad outputs in the
    // standard atmosphere calculation and plotting tabs

    computed = false;

    // resets all units
    QList<QComboBox*> comboBoxes = ui->tab->findChildren<QComboBox*>();
    for (QComboBox* comboBox : std::as_const(comboBoxes)) {
        comboBox->setCurrentIndex(0);  // Set to first item
    }

    // clears lineEdits
    QList<QLineEdit*> lineEdits = ui->tab->findChildren<QLineEdit*>();
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


void STDAtmo::on_Reset_graph_clicked()
{
    // clears the plot and removes saved data

    ui -> PlottingWidget -> clearGraphs();
    ui->PlottingWidget->xAxis->setLabel(" ");
    ui->PlottingWidget->yAxis->setLabel("Altitude, m");
    ui->PlottingWidget->xAxis->setRange(0,1);
    ui->PlottingWidget->yAxis->setRange(0,86000);
    ui->PlottingWidget->replot();

}

void STDAtmo::on_Reset_Yp_clicked()
{
    // removes all inputs nad outputs from the Y+ calculation list
    // and returns to default units

    // resets all units
    QList<QComboBox*> comboBoxes = ui->tab_2->findChildren<QComboBox*>();
    for (QComboBox* comboBox : std::as_const(comboBoxes)) {
        comboBox->setCurrentIndex(0);  // Set to first item
    }

    // clears lineEdits
    QList<QLineEdit*> lineEdits = ui->tab_2->findChildren<QLineEdit*>();
    for (QLineEdit* lineEdit : lineEdits) {
        lineEdit->clear();
    }

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

    // import class constants
    const auto& Href = standardAtmosphere::Href;


    // Read the input values and convert to double
    QString dISAtext = ui -> dISAInp -> text();
    QString Alttext = ui -> AltitudeInp -> text();

    // check if all inputs are written correctly
    if (dISAtext.isEmpty() || Alttext.isEmpty()) {
        QMessageBox::warning(this, "Error", "The input is empty. Enter the altitude and temperature deviation");
        return;
    }

    if (!isDouble(dISAtext) ) {
       QMessageBox::warning(this, "Error", "Enter numbers as inputs");
       ui->dISAInp->clear();
       return;
    }
    else if (!isDouble(Alttext)){
        QMessageBox::warning(this, "Error", "Enter numbers as inputs");
        ui->AltitudeInp->clear();
        return;
    }

    dISA = dISAtext.toDouble();
    H  = Alttext.toDouble();            // geometric altitude

    // compute a geopotential altitude
    Hgp = standardAtmosphere::compute_geopotential_altitude(H);

    // Read the units and convert to SI
    QString AltUnits = ui -> AltUnits -> currentText();
    QString dISAUnits = ui -> TempUnits -> currentText();

    if (AltUnits == "ft") { Hgp /= 3.28084; }
    if (dISAUnits == "°F / °R") { dISA = dISA * 5/9; }

    // Check the altitude limit

    if (Hgp < 0){
        QMessageBox::warning(this, "Error", "The altitude range must be greater than 0. Add a positive altitude");
        return;
    }

    if (Hgp > Href[7]){
        if (Hgp > Href[7]+ Heps0)     // adds a 0.5 m tolerance
        {
            QMessageBox::information(this,"Title","You are aiming too high, son...\n Enter an altitude below 86 km");
        }
        else{
            Hgp = Href[7];
        }
    }


    standardAtmosphere::compute_standard_atmosphere(Hgp, dISA, g, p, T, ro, mu, a);
    computed = true;

    // output results in appropriate units
    output_STDAresults();

}

void STDAtmo::on_GopAltCombo_currentIndexChanged(){output_STDAresults();}
void STDAtmo::on_PressureCombo_currentIndexChanged(){output_STDAresults();}
void STDAtmo::on_TemperatureCombo_currentIndexChanged(){output_STDAresults();}
void STDAtmo::on_DensityCombo_currentIndexChanged(){output_STDAresults();}
void STDAtmo::on_SpeedSoundCombo_currentIndexChanged(){output_STDAresults();}
void STDAtmo::on_ViscosityCombo_currentIndexChanged(){output_STDAresults();}
void STDAtmo::on_GravityCombo_currentIndexChanged(){output_STDAresults();}

void STDAtmo::output_STDAresults(){

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

    // Dynamic viscosity
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

void STDAtmo::convert_input_values_Ypp(double &Uinf, double &rhoinf, double &muinf, double &Lref){

    // converts all inputs to SI units depending on the local combobox units setup

    int selectedText;

    // Density
    selectedText = ui->DensityCombo_Yp->currentIndex();
    switch (selectedText) {
    case 0: break;
    case 1: rhoinf = Units::slugPerFt3ToKgPerM3(rhoinf); break;
    }

    // Airspeed
    selectedText = ui->AirspeedUnits_Yp->currentIndex();
    switch (selectedText) {
    case 0: break;
    case 1: Uinf = Units::kilometersPerHourToMetersPerSecond(Uinf); break;
    case 2: Uinf = Units::feetPerSecondToMetersPerSecond(Uinf); break;
    case 3: Uinf = Units::milesPerHourToMetersPerSecond(Uinf); break;
    case 4: Uinf = Units::KnotsToMetersPerSecond(Uinf); break;
    }

    // Dynamic viscosity
    selectedText = ui->ViscosityCombo_Yp->currentIndex();
    switch (selectedText) {
    case 0: break;
    case 1: muinf = Units::lbfSecondPerFt2ToPascalSecond(muinf); break;
    }

    // Characteristic length
    selectedText = ui->AltUnits_Yp->currentIndex();
    switch (selectedText) {
    case 0: break;
    case 1: muinf = Units::feetToMeter(muinf); break;
    }

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

    // import class constants
    const auto& Href = standardAtmosphere::Href;

    ui -> PlottingWidget -> clearGraphs();

    // plots standard atmosphere properties
    QString HminText = ui -> AltitudeMinInp -> text();
    QString HmaxText = ui -> AltitudeMaxInp -> text();
    QString dHText   = ui -> DAltInp -> text();
    QString pl_type  = ui -> UnitTypePlot -> currentText();
    QString unit_inp = ui -> UnitPlot -> currentText();
    QString unit_def = ui -> UnitSetup -> currentText();

    // check if any input is not a number
    if (!isDouble(HminText) ) {
        QMessageBox::warning(this, "Error", "Enter numbers as inputs");
        ui->AltitudeMinInp->clear();
        return;
    }
    else if (!isDouble(HmaxText)){
        QMessageBox::warning(this, "Error", "Enter numbers as inputs");
        ui->AltitudeMaxInp->clear();
        return;
        }
    else if (!isDouble(dHText)){
        QMessageBox::warning(this, "Error", "Enter numbers as inputs");
        ui->DAltInp->clear();
        return;
        }


    hbegin = HminText.toDouble();
    hend   = HmaxText.toDouble();
    dh     = dHText.toDouble();

    if (hbegin < 0 || hend > 86000 || dh < 0){
        QMessageBox::warning(this, "Error", "The altitude range must be between 0 and 86000 m with a "
                                            "positive altitude increment. Enter a valid region");
        return;
    }


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

    atmoData.clear();

    // compute geopotential altitudes
    double hbeg_geo = standardAtmosphere::compute_geopotential_altitude(hbegin);
    double hend_geo = standardAtmosphere::compute_geopotential_altitude(hend);

    int Npoints = static_cast<int>(std::ceil((hend - hbegin) / dh)) + 1;        // number of points to plot except for the lazers boundaries
    atmoData.reserve(Npoints + 8);   // reserves extra 8 points for potential layers

    // compute points to plot
    for (double X = hbegin; X <= hend + 1e-6; X += dh) {
        double Hgp = standardAtmosphere::compute_geopotential_altitude(X);
        standardAtmosphere::compute_standard_atmosphere(Hgp, dISA, g, p, T, ro, mu, a);
        append_selected_unit_value(standardAtmosphere::compute_geometric_altitude(Hgp), g, p, T, ro, mu, a, atmoData,
                                   pl_type, unit_inp, unit_def);
    }

    // determine how many layers and points are required
    for (int i = 0; i < 8; ++i) {
        if (Href[i] > hbeg_geo && Href[i] < hend_geo) {
            standardAtmosphere::compute_standard_atmosphere(Href[i], dISA, g, p, T, ro, mu, a);
            append_selected_unit_value(standardAtmosphere::compute_geometric_altitude(Href[i]), g, p, T, ro, mu, a, atmoData,
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
                    gp_alt = standardAtmosphere::compute_geopotential_altitude(Units::feetToMeter(atmoData[datacount].second));
                }
                else{
                    gp_alt = standardAtmosphere::compute_geopotential_altitude(atmoData[datacount].second);
                }
                // append calculated values into separate vectors
                xData.append(atmoData[datacount].second);
                yData.append(atmoData[datacount].first);
                datacount += 1;
            } while (gp_alt + Heps0 < Href[i] && datacount < data_size);         // used tolerance to avoid issue while converting units

            // plot the graph
            if (!xData.isEmpty()) {
                ui->PlottingWidget->addGraph();
                ui->PlottingWidget->graph(i-1)->setData(yData, xData);

                Xmax.append(*std::max_element(yData.begin(), yData.end()));
                Ymax.append(*std::max_element(xData.begin(), xData.end()));
                Xmin.append(*std::min_element(yData.begin(), yData.end()));
                Ymin.append(*std::min_element(xData.begin(), xData.end()));

                ui->PlottingWidget->graph()->setPen(QPen(Qt::blue, 2));
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


void STDAtmo::on_AltUnits_Yp_2_currentIndexChanged()
{
    // updates the first mesh layer step size units output
    if (!computed) return;

    int selectedText;

    // Altitude output
    selectedText = ui-> AltUnits_Yp_2 -> currentIndex();
    switch (selectedText){
    case 0: ui -> dsOutp_Yp -> setText(QString::number(ds, 'e', 4)); break;
    case 1: ui -> dsOutp_Yp -> setText(QString::number(Units::metersToFeet(ds), 'e', 4)); break;
    }
}


void STDAtmo::on_Compute_Yp_clicked()
{
    // computs the first mesh layer step size using a flat plate formulation

    // Extract inputs
    QString Vinf_text  = ui -> AirspeedInp_Yp -> text();
    QString roinf_text = ui -> DensityInp_Yp -> text();
    QString muinf_text = ui -> ViscosityInp_Yp -> text();
    QString Lref_text = ui -> LengthInp_Yp -> text();
    QString Yp_text = ui -> YpInp_Yp -> text();

    // check if all inputs are written correctly
    if (Vinf_text.isEmpty() || roinf_text.isEmpty() || muinf_text.isEmpty() || Lref_text.isEmpty() || Yp_text.isEmpty()) {
        QMessageBox::warning(this, "Error", "The input is empty. Enter all necessary inputs");
        return;
    }

    if (!isDouble(Vinf_text) ) {
        QMessageBox::warning(this, "Error", "Enter numbers as inputs");
        ui->AirspeedInp_Yp->clear();
        return;
    }
    else if (!isDouble(roinf_text)){
        QMessageBox::warning(this, "Error", "Enter numbers as inputs");
        ui->DensityInp_Yp->clear();
        return;
    }
    else if (!isDouble(muinf_text)){
        QMessageBox::warning(this, "Error", "Enter numbers as inputs");
        ui->ViscosityInp_Yp->clear();
        return;
    }
    else if (!isDouble(Lref_text)){
        QMessageBox::warning(this, "Error", "Enter numbers as inputs");
        ui->LengthInp_Yp->clear();
        return;
    }
    else if (!isDouble(Yp_text)){
        QMessageBox::warning(this, "Error", "Enter numbers as inputs");
        ui->YpInp_Yp->clear();
        return;
    }


    // convert inputs into SI units
    double rho_inf = roinf_text.toDouble();
    double Uinf = Vinf_text.toDouble();
    double muinf = muinf_text.toDouble();
    double Lref = Lref_text.toDouble();
    double Ypp  = Yp_text.toDouble();

    if (rho_inf < 0 || Uinf < 0 || muinf < 0 || Lref < 0 || Ypp < 0) {
        QMessageBox::warning(this, "Error", "All inputs must be positive");
        return;
    }

    convert_input_values_Ypp(Uinf, rho_inf, muinf, Lref);

    // compute Re and first step size
    double Re    = rho_inf * Uinf * Lref / muinf;
    double Cf    = 0.455 / pow(log10(Re), 2.58);
    double tau_w = 0.5 * Cf * rho_inf * pow(Uinf, 2);
    double Ufric = pow(tau_w/rho_inf, 0.5);
    ds = Ypp * muinf / (Ufric * rho_inf);

    // output results
    ui -> ReOutp_Yp -> setText(QString::number(Re, 'e', 3));

    QString output_units = ui -> AltUnits_Yp_2 ->currentText();
    if (output_units == "ft"){
        ui -> dsOutp_Yp -> setText(QString::number(Units::metersToFeet(ds), 'e', 3));
    }
    else{
        ui -> dsOutp_Yp -> setText(QString::number(ds, 'e', 4));
    }

}


void STDAtmo::on_Reset_SP_clicked()
{
    // resets all units and clears all inputs
    // and outputs in the airspeed tab

    computed = false;

    // resets all units
    QList<QComboBox*> comboBoxes = ui->tab_5->findChildren<QComboBox*>();
    for (QComboBox* comboBox : std::as_const(comboBoxes)) {
        comboBox->setCurrentIndex(0);  // Set to first item
    }

    // clears lineEdits
    QList<QLineEdit*> lineEdits = ui->tab_5->findChildren<QLineEdit*>();
    for (QLineEdit* lineEdit : lineEdits) {
        lineEdit->clear();
    }


}


void STDAtmo::on_Compute_SP_clicked()
{

    // computes airspeeds based on the given speed,
    // altitude and temperature devaition inputs

    // import class constants
    const auto& Href = standardAtmosphere::Href;


    // Read the input values and convert to double
    QString dISAtext = ui -> dISAInp_SP -> text();
    QString Alttext = ui -> AltitudeInp_SP -> text();

    // check if all inputs are written correctly
    if (dISAtext.isEmpty() || Alttext.isEmpty()) {
        QMessageBox::warning(this, "Error", "The input is empty. Enter the values");
        return;
    }

    if (!isDouble(dISAtext) ) {
        QMessageBox::warning(this, "Error", "Enter numbers as inputs");
        ui->dISAInp_SP->clear();
        return;
    }
    else if (!isDouble(Alttext)){
        QMessageBox::warning(this, "Error", "Enter numbers as inputs");
        ui->AltitudeInp_SP->clear();
        return;
    }


    dISA = dISAtext.toDouble();
    H  = Alttext.toDouble();            // geometric altitude

    if (H < 0 ) {
        QMessageBox::warning(this, "Error", "Enter a positive altitude");
        return;
    }

    // compute a geopotential altitude
    Hgp = standardAtmosphere::compute_geopotential_altitude(H);

    // Read the units and convert to SI
    QString AltUnits = ui -> AltUnits_SP -> currentText();
    QString dISAUnits = ui -> TempUnits_SP -> currentText();

    if (AltUnits == "ft") { Hgp /= 3.28084; }
    if (dISAUnits == "°F / °R") { dISA = dISA * 5/9; }

    // Check the altitude limit
    if (Hgp > Href[7]){
        if (Hgp > Href[7]+ Heps0)     // adds a 0.5 m tolerance
        {
            QMessageBox::information(this,"Title","You are aiming too high, son...\n Enter an altitude below 86 km");
        }
        else{
            Hgp = Href[7];
        }
    }

    standardAtmosphere::compute_standard_atmosphere(Hgp, dISA, g, p, T, ro, mu, a);


    // convert airspeeds
    bool CAS_flag = ui -> CASradioButton->isChecked();
    bool EAS_flag = ui -> EASradioButton->isChecked();
    bool TAS_flag = ui -> TASradioButton->isChecked();
    bool Mach_flag = ui -> MachradioButton->isChecked();

    if (TAS_flag == 1){


        // extract TAS
        QString VTAS_text = ui -> TASInpOut -> text();

        // check if the airspeed was entered
        if (VTAS_text.isEmpty()) {
            QMessageBox::warning(this, "Error", "The input is empty. Enter the airspeed");
            return;
        }

        // convert to m/s
        VTAS = VTAS_text.toDouble();
        checkAirspeed(VTAS);
        VTAS = STDAtmo::convert_airspeed_Input(ui -> TASCombo, VTAS);

        // compute EAS
        VEAS = standardAtmosphere::EASfromTAS(ro, VTAS);

        // compute Mach number
        Mach = standardAtmosphere::MachfromTAS(a, VTAS);

        // compute calibrated airspeed
        VCAS = standardAtmosphere::CASfromEAS(p, VEAS, Mach);

    }

    else if (CAS_flag == 1){

        // extract CAS
        QString VCAS_text = ui -> CASInpOut -> text();

        // check if the airspeed was entered
        if (VCAS_text.isEmpty()) {
            QMessageBox::warning(this, "Error", "The input is empty. Enter the airspeed");
            return;
        }

        // convert to m/s
        VCAS = VCAS_text.toDouble();
        checkAirspeed(VCAS);
        VCAS = STDAtmo::convert_airspeed_Input(ui -> CASCombo, VCAS);

        // compute Mach number
        Mach = standardAtmosphere::MachfromCAS(p, VCAS);

        // compute TAS
        VTAS = standardAtmosphere::TASfromMach(a, Mach);

        // compute EAS
        VEAS = standardAtmosphere::EASfromTAS(ro, VTAS);

    }

    else if (EAS_flag == 1){

        // extract EAS
        QString VEAS_text = ui -> EASInpOut -> text();

        // check if the airspeed was entered
        if (VEAS_text.isEmpty()) {
            QMessageBox::warning(this, "Error", "The input is empty. Enter the airspeed");
            return;
        }

        // convert to m/s
        VEAS = VEAS_text.toDouble();
        checkAirspeed(VEAS);
        VEAS = STDAtmo::convert_airspeed_Input(ui -> EASCombo, VEAS);

        // compute TAS
        VTAS = standardAtmosphere::TASfromEAS(ro, VEAS);

        // compute Mach number
        Mach = standardAtmosphere::MachfromTAS(a, VTAS);

        // compute CAS
        VCAS = standardAtmosphere::CASfromEAS(p, VEAS, Mach);

    }

    else if (Mach_flag == 1){

        // extract Mach
        QString Mach_text = ui -> MachInpOut -> text();
        checkAirspeed(Mach);
        Mach = Mach_text.toDouble();

        // check if the airspeed was entered
        if (Mach_text.isEmpty()) {
            QMessageBox::warning(this, "Error", "The input is empty. Enter the Mach number");
            return;
        }

        // compute TAS
        VTAS = standardAtmosphere::TASfromMach(a, Mach);

        // compute EAS
        VEAS = standardAtmosphere::EASfromTAS(ro, VTAS);

        // compute CAS
        VCAS = standardAtmosphere::CASfromEAS(p, VEAS, Mach);

    }

    computed = true;

    // output results
    convert_airspeed(ui->CASCombo, ui->CASInpOut, VCAS);
    convert_airspeed(ui->TASCombo, ui->TASInpOut, VTAS);
    convert_airspeed(ui->EASCombo, ui->EASInpOut, VEAS);

    ui -> MachInpOut -> setText(QString::number(Mach, 'g', 3));

}


void STDAtmo::on_CASradioButton_clicked()
{
    // adjust lineEdits. CAS
    ui->CASInpOut->setReadOnly(false);
    ui->CASInpOut->setStyleSheet("");
    ui->TASInpOut->setReadOnly(true);
    ui->TASInpOut->setStyleSheet("background-color: #222222;");
    ui->EASInpOut->setReadOnly(true);
    ui->EASInpOut->setStyleSheet("background-color: #222222;");
    ui->MachInpOut->setReadOnly(true);
    ui->MachInpOut->setStyleSheet("background-color: #222222;");
}


void STDAtmo::on_MachradioButton_clicked()
{
    // adjust lineEdits. Mach
    ui->TASInpOut->setReadOnly(true);
    ui->TASInpOut->setStyleSheet("background-color: #222222;");
    ui->EASInpOut->setReadOnly(true);
    ui->EASInpOut->setStyleSheet("background-color: #222222;");
    ui->CASInpOut->setReadOnly(true);
    ui->CASInpOut->setStyleSheet("background-color: #222222;");
    ui->MachInpOut->setReadOnly(false);
    ui->MachInpOut->setStyleSheet("");
}


void STDAtmo::on_EASradioButton_clicked()
{
    // adjust lineEdits. EAS
    ui->TASInpOut->setReadOnly(true);
    ui->TASInpOut->setStyleSheet("background-color: #222222;");
    ui->MachInpOut->setReadOnly(true);
    ui->MachInpOut->setStyleSheet("background-color: #222222;");
    ui->CASInpOut->setReadOnly(true);
    ui->CASInpOut->setStyleSheet("background-color: #222222;");
    ui->EASInpOut->setReadOnly(false);
    ui->EASInpOut->setStyleSheet("");
}


void STDAtmo::on_TASradioButton_clicked()
{
    // adjust lineEdits. TAS
    ui->EASInpOut->setReadOnly(true);
    ui->EASInpOut->setStyleSheet("background-color: #222222;");
    ui->MachInpOut->setReadOnly(true);
    ui->MachInpOut->setStyleSheet("background-color: #222222;");
    ui->CASInpOut->setReadOnly(true);
    ui->CASInpOut->setStyleSheet("background-color: #222222;");
    ui->TASInpOut->setReadOnly(false);
    ui->TASInpOut->setStyleSheet("");
}


void STDAtmo::on_TASCombo_currentIndexChanged(){convert_airspeed(ui-> TASCombo, ui-> TASInpOut, VTAS);}
void STDAtmo::on_EASCombo_currentIndexChanged(){convert_airspeed(ui-> EASCombo, ui-> EASInpOut, VEAS);}
void STDAtmo::on_CASCombo_currentIndexChanged(){convert_airspeed(ui-> CASCombo, ui-> CASInpOut, VCAS);}


void STDAtmo::on_Export_graph_clicked()
{
    // export plotted results into a CSV file
    std::ofstream outputFile;

    QString fileName = QFileDialog::getSaveFileName(this,
                                                    "Save Plot Data",                      // Dialog title
                                                    QDir::homePath() + "/STDAtmo_plot.csv", // Default directory and filename
                                                    "CSV Files (*.csv);;All Files (*)"      // File filters
                                                    );

    if (!fileName.isEmpty()) {
        QFile file(fileName);
        if (file.open(QIODevice::WriteOnly | QIODevice::Text)) {
            // Use QTextStream instead of std::ofstream
            QTextStream out(&file);

            // Write header - convert QString to stream
            out << "Altitude, " << plot_unitsY << " ; "
                << plot_unit_type << ", " << plot_unitsX << "\n";

            // Write data
            for (int i = 0; i < atmoData.size(); ++i) {
                out << atmoData[i].second << " ; " << atmoData[i].first << "\n";
            }

            file.close();

            QMessageBox::information(this, "Success",
                                     QString("File saved to:\n%1").arg(fileName));
        } else {
            QMessageBox::warning(this, "Error",
                                 QString("Could not save file:\n%1\nError: %2")
                                     .arg(fileName, file.errorString()));
        }
    }

}

