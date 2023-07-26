#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <string>
#include <sys/stat.h>
#include <time.h>
#include "stdio.h"
#include <unistd.h>
#include <cstdlib>

#include <TCanvas.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TH2F.h>
#include <TF1.h>
#include <TH1F.h>
#include <TRandom.h>
#include <TMath.h>

#include "Garfield/TrackHeed.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/SolidTube.hh"
#include "Garfield/GeometrySimple.hh"
#include "Garfield/ViewGeometry.hh"
#include "Garfield/ComponentConstant.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/FundamentalConstants.hh"
#include "Garfield/Random.hh"
#include "Garfield/Plotting.hh"
#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/AvalancheMC.hh"

using namespace Garfield;
using namespace std;

TRandom degrad_random;

/**
 * Execute a shell command
 */
string exec(const char* cmd){
    char buffer[128];
    string result = "";
    FILE* pipe = popen(cmd, "r");
    if (!pipe) throw runtime_error("popen() failed!");
    try {
        while (fgets(buffer, sizeof buffer, pipe) != NULL) {
            result += buffer;
        }
    } catch (...) {
        pclose(pipe);
        throw;
    }
    pclose(pipe);
    return result;
}

string fixedLength(int value){
    int digits = 6;
    string result;
    while (digits-- > 0){
        result += ('0' + value % 10);
        value /= 10;
    }
    if (value < 0){
        result += '-';
    }
    reverse(result.begin(), result.end());
    return result;
}

/**
 *  Returns a string containing the time information
 */
string time(){
    time_t Zeitstempel;
    tm *nun;
    Zeitstempel = time(0);
    nun = localtime(&Zeitstempel);
    int start = (nun->tm_hour) * 10000 + (nun->tm_min) * 100 + nun->tm_sec;
    string time = fixedLength(start);
    time.insert(2, "-");
    time.insert(5, "-");
    return time;
}

/**
 *  Returns a string containing the date in reverse order
 */
string date(){
    time_t Zeitstempel;
    tm *nun;
    Zeitstempel = time(0);
    nun = localtime(&Zeitstempel);
    int date = ((nun->tm_year + 1900) % 100) * 10000 + (nun->tm_mon + 1) * 100 + (nun->tm_mday);
    string date_s = fixedLength(date);
    return date_s;
}

/**
 *  Writes the run.txt file, To be called after (!!!) all the simulations
 */
void runfile(string begindate, string begintime, int runnumber, string dir){
    fstream f;
    string file = dir + "//run.txt";
    f.open(file, ios::out);
    f << "Run number run:\t\t" << runnumber << "\n";
    f << "Run start time:\t\t" << begindate << "_" << begintime << "\n";
    f << "Run end time:\t" << date() << "_" << time() << "\n";
    f << "\n";
    f << "Run parameters:\n";
    f << "\t run mode (0=time,1=number of frames) \t1\n";
    f << "\t run time / number of triggers      0\n";
    f << "\t shutter mode (0 = untiggered, 1 = untriggered long, 2 = untriggered very long, 3 = untiggered 2x faster clock, 4 = untiggered long 2x faster clock, 5 = untiggered very long 2x faster clock, 6 = external trigger, 7 = external trigger 2x faster clock   long\n";
    f << "\t shutter time 200\n";
    f << "\t readout mode (0 = zero suppressed, 1 = complete matrix)    0\n\n";
    f << "Chip IDs:\n";
    f << "\t FEC 0 Board 0 Chip 1: H1-W10\n";
}

/**
 *  Copies the config files provided in an extra folder into the run folder
 */
void makeconfig(string dir){
    //threshold
    string file = dir + "//chip_1_board_0_fec_0_threshold.txt";
    ifstream srce("config_vorlagen//chip_1_board_0_fec_0_threshold.txt", ios::binary);
    ofstream dest(file, ios::binary);
    dest << srce.rdbuf();
    //matrix
    file = dir + "//chip_1_board_0_fec_0_matrix.txt";
    srce = ifstream("config_vorlagen//chip_1_board_0_fec_0_matrix.txt", ios::binary);
    dest = ofstream(file, ios::binary);
    dest << srce.rdbuf();
    //fsr
    file = dir + "//chip_1_board_0_fec_0_fsr.txt";
    srce = ifstream("config_vorlagen//chip_1_board_0_fec_0_fsr.txt", ios::binary);
    dest = ofstream(file, ios::binary);
    dest << srce.rdbuf();
}

void write_text_to_position_file(const char* file, const char* text){
    FILE* pFile = fopen(file, "a");
    fprintf(pFile,"%s\n", text);
    fclose(pFile);
}

void generateGasFile(std::string filename, std::string gas1, double fraction1, std::string gas2, double fraction2, double pressure, double temperature, double ncollisions, double emin, double emax, unsigned int nE, double bmin, double bmax, unsigned int nB) {
  // Define the gas mixture
  Garfield::MediumMagboltz gas(gas1.c_str(), fraction1, gas2.c_str(), fraction2);
  gas.SetTemperature(temperature);
  gas.SetPressure(pressure);
  // Define the field grid
  gas.SetFieldGrid(emin, emax, nE, false, bmin, bmax, nB);
  // Generate the gas table
  gas.GenerateGasTable(ncollisions);
  // Generate the gas file
  gas.WriteGasFile(filename.c_str());
}

/**
 *  Simulate the conversion point of a photon with a given energy and a start point
 *  Use a defined Garfield sensor (with gas mixture, temperature and pressure)
 */
double get_start_position(double energy, Sensor& sensor, double z_start){
    // Initialize the photon track
    TrackHeed trackphoton;
    trackphoton.SetSensor(&sensor);
    trackphoton.EnableElectricField();

    // Initial coordinates of the photon.
    const double x0 = 0.;
    const double y0 = 0.;
    const double z0 = z_start;
    const double t0 = 0.;
    int ne = 0;

    // Simulate the photon track
    trackphoton.TransportPhoton(x0, y0, z0, t0, energy, 0., 0., -1., ne);
    double z_sum = 0;

    // Iterate over the generated primary electrons to get their starting points
    for (int i = 0; i < ne; i++){
        // Variables for the start values of the electrons
        double start_x = 0.;
        double start_y = 0.;
        double start_z = 0.;
        double start_time = 0.;
        // Further variables - needed but not meaningful for this application
        double start_energy = 0.;
        double start_dx = 0.;
        double start_dy = 0.;
        double start_dz = 0.;

        // Fill the variables for the individual electrons
        trackphoton.GetElectron(i, start_x, start_y, start_z, start_time, start_energy, start_dx, start_dy, start_dz);
        z_sum += start_z;
    }
    if(ne == 0){
        return false;
    }

    // Get and return the mean starting position of the electrons
    double z_mean = z_sum / ne;
    return z_mean;
}

/**
 *  Simulate conversion points of photons with a given energy and a start point
 *  Use a defined Garfield sensor (with gas mixture, temperature and pressure)
 *  The function return 3 fit parameters of the resulting exponential curve
 */
void get_absoption_curve(string filename, double energy, Sensor& sensor, double z_start, double *fit_parameters){
    // Initialize the photon track
    TrackHeed trackphoton;
    trackphoton.SetSensor(&sensor);
    trackphoton.EnableElectricField();

    // Initial coordinates of the photon.
    const double x0 = 0.;
    const double y0 = 0.;
    const double z0 = z_start;
    const double t0 = 0.;
    int ne = 0;
    int total_photons = 0;

    // Create the histogram for the simulated absorption points
    TH1F absorption_points = TH1F("absorption_points", "absorption_points", 50, 0, z0);

    // Simulate 100000 photons
    for(int photon = 0; photon < 100000; photon++){
        // Simulate the photon track
        trackphoton.TransportPhoton(x0, y0, z0, t0, energy, 0., 0., -1., ne);
        double z_sum = 0;
        total_photons++;

        // If the photon is not absorbed within the gas volume simulate again
        if(ne == 0){
            photon--;
            continue;
        }
        // Iterate over the generated primary electrons to get their starting points
        for(int i = 0; i < ne; i++){
            // Variables for the start values of the electrons
            double start_x = 0.;
            double start_y = 0.;
            double start_z = 0.;
            double start_time = 0.;
            // Further variables - needed but not meaningful for this application
            double start_energy = 0.;
            double start_dx = 0.;
            double start_dy = 0.;
            double start_dz = 0.;

            // Fill the variables for the individual electrons
            trackphoton.GetElectron(i, start_x, start_y, start_z, start_time, start_energy, start_dx, start_dy, start_dz);
            z_sum += start_z;
        }

        // Get the mean starting position of the electrons and fill it into the histogram
        double z_mean = z_sum / ne;
        absorption_points.Fill(z_mean);
    }

    // Draw the histogram in a canvas
    auto c1 = new TCanvas();
    absorption_points.Draw();

    gStyle->SetOptStat(1);
    gStyle->SetOptFit(1);

    // Fit the histogram
    TF1 exp = TF1("fit","[0]*exp([1]*(x-[2]))", 0.3, 0.8*z0);
    exp.SetParLimits(1, 0, 10000);
    absorption_points.Fit("fit", "R");
    string title = "Absorption curve - Efficiency = " + to_string(100*100000.0/total_photons) + " %";
    absorption_points.SetTitle(title.c_str());
    absorption_points.Draw();

    // Return the fit parameters via a pointer
    fit_parameters[0] = exp.GetParameter(0);
    fit_parameters[1] = exp.GetParameter(1);
    fit_parameters[2] = exp.GetParameter(2);

    // Print the resulting histogram as PDF
    filename = filename + "_absorption.pdf";
    c1->Print(filename.c_str());

    return;
}

/**
 * Write an input file for degrad based on photon energy, the gas mixture (two gases and their percentages),
 * the gas pressure and temperature, the electric field.
 * It is written into a file with a given filename
 */
int write_degrad_file(double energy, int gas1, int gas2, double percentage1, double percentage2, double temperature, double pressure, double efield, string filename){
    // Get a random seed for degrad
    int seed = degrad_random.Integer(1000000);
    //cout << seed << endl;

    // As degrad needs the right amount of spaces fill the file with spaces based on the length of the numbers
    double parameters[9] = {(double)seed, energy, (double)gas1, (double)gas2, percentage1, percentage2, temperature, pressure, efield};
    int spaces[9] = {10, 6, 5, 5, 6, 6, 7, 7, 7};
    string space_strings[9] = {" ", " ", " ", " ", " ", " ", " ", " ", " "};

    for(int i = 0; i < 9; i++){
        int comparison = 1;
        while(true){
            if(parameters[i] >= comparison){
                spaces[i]--;
                comparison *= 10;
            }
            else{
                break;
            }
        }
        for(int j = 1; j < spaces[i]; j++){
            space_strings[i] += " ";
        }
    }

    // Open and fill the file
    ofstream infile;
    infile.open(filename);
    infile << "         2        10         3        -1" << fixed << setprecision(0) << space_strings[0] << parameters[0] << setprecision(3) << space_strings[1] << parameters[1] << "     2.0000   100.0000\n";
    infile << fixed << setprecision(0) << space_strings[2] << parameters[2] << space_strings[3] << parameters[3] << "   99   99   99   99\n";
    infile << fixed << setprecision(3) << space_strings[4] << parameters[4] << space_strings[5] << parameters[5] << "     0.000     0.000     0.000     0.000" << space_strings[6] << parameters[6] << space_strings[7] << parameters[7] << "\n";
    infile << fixed << setprecision(2) << space_strings[8] << parameters[8] << "      0.00      0.00    1    1\n";
    infile << "   100.00      0.500    0    0    1    1    1    1    1\n";
    infile << "\n";
    infile << "\n";
    infile.close();

    return seed;
}

/**
 * Run degrad with a given input file
 */
void run_degrad(string filename){
    string cmd_rm;
    string cmd_rm2;
    string cmd_degrad;
    cmd_rm = "rm -f DEGRAD.OUT";
    cmd_rm2 = "rm -f degrad.out";
    cmd_degrad = "./degrad.run < " + filename + " > degrad.out";
    //cout << cmd_rm << "\t" << cmd_rm2 << "\t" << cmd_degrad << endl;
    string err = "";

    // Remove potential *.OUT files
    err = exec(cmd_rm.c_str());
    //cout << cmd_rm << ":\t" << err << endl; // Debug output

    // Remove potential *.out files
    err = exec(cmd_rm2.c_str());
    //cout << cmd_rm2 << ":\t" << err << endl; // Debug output

    // Start degrad
    err = exec(cmd_degrad.c_str());
    //cout << cmd_degrad << ":\t" << err << endl; // Debug output
}

/**
 * Get the degrad integer for a gas based on a Garfield++ gas name
 */
int get_gas_parameter(string gasname){
    if (gasname == "CF4"){return 1;}
    else if (gasname == "Ar"){return 2;}
    else if (gasname == "Argon"){return 2;}
    else if (gasname == "He"){return 3;}
    else if (gasname == "Helium"){return 3;}
    else if (gasname == "He3"){return 4;}
    else if (gasname == "Helium-3"){return 4;}
    else if (gasname == "Ne"){return 5;}
    else if (gasname == "Neon"){return 5;}
    else if (gasname == "Kr"){return 6;}
    else if (gasname == "Krypton"){return 6;}
    else if (gasname == "Xe"){return 7;}
    else if (gasname == "Xenon"){return 7;}
    else if (gasname == "CH4"){return 8;}
    else if (gasname == "Methane"){return 8;}
    else if (gasname == "Ethane"){return 9;}
    else if (gasname == "Propane"){return 10;}
    else if (gasname == "Isobutane"){return 11;}
    else if (gasname == "CO2"){return 12;}
    else if (gasname == "Neo-Pentane"){return 13;}
    else if (gasname == "H2O"){return 14;}
    else if (gasname == "Water"){return 14;}
    else if (gasname == "O"){return 15;}
    else if (gasname == "Oxygen"){return 15;}
    else if (gasname == "N"){return 16;}
    else if (gasname == "Nitrogen"){return 16;}
    else if (gasname == "Nitric Oxide"){return 17;}
    else if (gasname == "Nitrous Oxide"){return 18;}
    else if (gasname == "Ethene"){return 19;}
    else if (gasname == "Acetylene"){return 20;}
    else if (gasname == "H2"){return 21;}
    else if (gasname == "Hydrogen"){return 21;}
    else if (gasname == "D2"){return 22;}
    else if (gasname == "Deuterium"){return 22;}
    else if (gasname == "CO"){return 23;}
    else if (gasname == "Carbon Monoxide"){return 23;}
    else if (gasname == "Methylal"){return 24;}
    else if (gasname == "DME"){return 25;}
    else if (gasname == "Reid Step Model"){return 26;}
    else if (gasname == "Maxwell Model"){return 27;}
    else if (gasname == "Reid Ramp Model"){return 28;}
    else if (gasname == "C2F6"){return 29;}
    else if (gasname == "SF6"){return 30;}
    else if (gasname == "NH3"){return 31;}
    else if (gasname == "Ammonia"){return 31;}
    else if (gasname == "C3H6"){return 32;}
    else if (gasname == "Propene"){return 32;}
    else if (gasname == "Cylopropane"){return 33;}
    else if (gasname == "CH3OH"){return 34;}
    else if (gasname == "Methanol"){return 34;}
    else if (gasname == "C2H5OH"){return 35;}
    else if (gasname == "Ethanol"){return 35;}
    else if (gasname == "C3H7OH"){return 36;}
    else if (gasname == "Iso-Propanol"){return 36;}
    else if (gasname == "Cs"){return 37;}
    else if (gasname == "Cesium"){return 37;}
    else if (gasname == "F"){return 38;}
    else if (gasname == "Flourine"){return 38;}
    else if (gasname == "CS2"){return 39;}
    else if (gasname == "COS"){return 40;}
    else if (gasname == "CD4"){return 41;}
    else if (gasname == "BF3"){return 42;}
    else if (gasname == "Boron-Triflouride"){return 42;}
    else if (gasname == "C2HF5"){return 43;}
    else if (gasname == "C2H2F4"){return 43;}
    else if (gasname == "TMA"){return 44;}
    else if (gasname == "C3H7OH"){return 46;}
    else if (gasname == "N-Propanol"){return 46;}
    else if (gasname == "CHF3"){return 50;}
    else if (gasname == "CF3BR"){return 51;}
    else if (gasname == "C3F8"){return 52;}
    else if (gasname == "O3"){return 53;}
    else if (gasname == "Ozone"){return 53;}
    else if (gasname == "Hg"){return 54;}
    else if (gasname == "Mercury"){return 54;}
    else if (gasname == "H2S"){return 55;}
    else if (gasname == "N-Butane"){return 56;}
    else if (gasname == "N-Pentane"){return 57;}
    else if (gasname == "CCL4"){return 61;}
    else {return 1;}
}

int main(int argc, char *argv[]){
    // Get the timestamp - needed for folder- and filenames
    string dateb = date();
    string timeb = time();

    // Initialize the pixelmatrix
    double pixelsize = 0.0055;
    const int pixel = 256; //in one direction
    int hits[pixel][pixel];

    bool create_gasfile = true;
    string gasfile;

    if (argc < 20){
        cout << "There are missing some arguments. The command is ./simulation <path> <job> <absorption> <approach> <length> <energy> <gas1> <gas2> <percentage1> <percentage2> <temperature> <pressure> <field> <polarization> <angle_offset> <amp_scaling> <amp_gain> <amp_width> <events> <degrad_output>" << endl;
        cout << "If no gasfile is provided a new one is generated (takes a couple of hours)" << endl;
        return 1;
    }
    if (argc == 21){
        gasfile = argv[1];
        create_gasfile = false;
    }
    else{
        create_gasfile = true;
    }
    int job = atoi(argv[argc - 19]);
    int absorption_approach = atoi(argv[argc - 18]);
    int simulation_approach = atoi(argv[argc - 17]);
    double length = atof(argv[argc - 16]);
    double energy = atof(argv[argc - 15]);
    string gas1 = argv[argc - 14];
    string gas2 = argv[argc - 13];
    double percentage1 = atof(argv[argc - 12]);
    double percentage2 = atof(argv[argc - 11]);
    double temperature = atof(argv[argc - 10]);
    double pressure = atof(argv[argc - 9]);
    double efield = atof(argv[argc - 8]);
    double polarization = atof(argv[argc - 7]);
    double angle_offset = atof(argv[argc - 6]);
    double amplification_scaling = atof(argv[argc - 5]);
    double amplification_gain = atof(argv[argc - 4]);
    double amplification_width = atof(argv[argc - 3]);
    int nEvents = atoi(argv[argc - 2]);
    int degrad_output = atoi(argv[argc - 1]);

    if (create_gasfile){
        cout << "MAGBOLTZ: Generate gasfile" << endl;
        gasfile = gas1 + "_" + gas2 + "_" + to_string(int(percentage1)) + "_" + to_string(int(percentage2)) + "_" + to_string(int(pressure)) + "torr_" + to_string(int(temperature)) + "C_" + to_string(int(efield)) + "Vcm.gas";
        generateGasFile(gasfile, gas1, percentage1, gas2, percentage2, pressure, 273.15+temperature, 20, efield, efield, 1, 0,0,1);
    }

    // Create the needed folders
    string dir = "Run_" + fixedLength(job) + "_" + date() + "_" + time();
    mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    makeconfig(dir);
    string degrad_dir = "Run_" + fixedLength(job) + "_" + date() + "_" + time() + "_degrad";
    if(degrad_output == 1){
        mkdir(degrad_dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    }

    // Make a gas medium.
    std::unique_ptr<MediumMagboltz> gas = std::make_unique<MediumMagboltz>();
    gas->SetComposition(gas1, percentage1, gas2, percentage2);
    gas->SetTemperature(273.15+temperature);
    gas->SetPressure(pressure);
    gas->LoadGasFile(gasfile);
    gas->Initialise(true);

    // Get diffusion parameters and drift velocity
    double el_Tdiff = 0.0;
    gas->GetElectronTransverseDiffusion(0, 0, 0, el_Tdiff);
    double el_Ldiff = 0.0;
    gas->GetElectronLongitudinalDiffusion(0, 0, 0, el_Ldiff);
    double el_vel = 0.0;
    gas->GetElectronVelocityE(0, 0, 0, el_vel);

    // Create a cylinder in which the x-rays can convert.
    // Diameter [cm]
    const double diameter = 7.8;
    SolidTube tube(0.0, 0.0, length, 0.5 * diameter, length);

    // Combine gas and box to a simple geometry.
    GeometrySimple geo;
    geo.AddSolid(&tube, gas.get());

    //Debug plots:
    //ViewGeometry geoView;
    //geoView.SetGeometry(&geo);
    //geoView.Plot3d();
    //geoView.SetPlaneXY();
    //geoView.Plot2d();
    //app.Run(true);

    // Make a component with constant electric field.
    ComponentConstant field;
    field.SetGeometry(&geo);
    field.SetElectricField(0., 0., efield);

    // Make a sensor.
    Sensor sensor;
    sensor.AddComponent(&field);

    TF1 absoption_curve = TF1("absoption_curve","[0]*exp([1]*(x-[2]))", 0, length);

    // If the absorption is simulated via monte carlo get the absorption curve via Garfield
    if(absorption_approach != 0){
        cout << "GARFIELD: Simulate 100000 photons to get the absorption curve" << endl;
        double fit_parameters[3];
        get_absoption_curve(dir, energy, std::ref(sensor), length, fit_parameters);
        absoption_curve.SetParameter(0, fit_parameters[0]);
        absoption_curve.SetParameter(1, fit_parameters[1]);
        absoption_curve.SetParameter(2, fit_parameters[2]);
        cout << "GARFIELD: Absorption curve simulation finished" << endl;
    }

    int event = 0;
    int photoelectrons = 0;
    while (true){
        int number = 0; //of activated pixels

        // Clean up the pixelmatrix for a new event
        for (int l = 0; l < pixel; l++){
            for (int m = 0; m < pixel; m++){
                hits[l][m] = 0;
            }
        }

        cout << "GARIFELD: photon " << event << endl;

        TF1 angles = TF1("angles","[0]*TMath::Cos(x) * TMath::Cos(x)+[1]", 0, 2*TMath::Pi());
        angles.SetParameter(0, 1 - (1 - polarization) / (1 + polarization));
        angles.SetParameter(1, (1 - polarization) / (1 + polarization));
        double angle = angles.GetRandom() + angle_offset;

        // Generate the filenames for the current event
        string filename = "/run_" + fixedLength(job) + "_data_" + fixedLength(event) + "_" + date() + "_" + time();
        string in_file;
        if(degrad_output == 1){
            in_file = degrad_dir + filename + ".in";
        }
        else{
            in_file = "degrad.in";
        }

        // Create the degrad in file and run degrad with it
        int degrad_seed;
        degrad_seed = write_degrad_file(energy, get_gas_parameter(gas1), get_gas_parameter(gas2), percentage1, percentage2, temperature, pressure, efield, in_file);
        run_degrad(in_file);
        if(degrad_output != 1){
            string cmd_rm_in;
            cmd_rm_in = "rm -f degrad.in";
            string err = "";
            err = exec(cmd_rm_in.c_str());
        }
        cout << "DEGRAD: Finished photoelectron track" << endl;

        double position;
        // Simulation of the absorption point per event via Garfield
        if(absorption_approach == 0){
            // Get a photon conversion point for the event. As an electron might not convert within the detector repeat until a position within the detector is found
            while (true){
                position = get_start_position(energy, std::ref(sensor), length);
                //cout << position << endl;
                if(position != false && position <= length && position >= 0.0){
                    break;
                }
            }
        }
        // Simulation of the absorption point via monte carlo and the absorption curve
        else{
            position = absoption_curve.GetRandom();
        }

        string file = dir + filename + ".txt";
        fstream f;
        f.open(file, ios::out);

        cout << "GARFIELD: Final position: " << position << endl;

        // Initialize variables for reading degrad data and for the secondary electrons
        Double_t x, y, z, t;
        Int_t n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n11;
        Int_t nlines = 0;
        ifstream in;
        Int_t nevt, nclus, nstexc, mcomp, mpair;
        double x1, y1, z1, t1, e1;
        double x2, y2, z2, t2, e2;
        int status;
        mpair = 0;

        // Open the degrad output file
        in.open("DEGRAD.OUT");

        //ViewDrift view; // For debugging

        // Read the event parameters from the degrad output file
        in >> nevt >> nclus >> nstexc >> mcomp >> n4 >> n5 >> n6 >> n7 >> n8 >> n9 >> n10 >> n11;
        printf("DEGRAD: nevt   = %d \n", nevt);
        printf("DEGRAD: nclus  = %d \n", nclus);
        printf("DEGRAD: nstexc = %d \n", nstexc);
        //printf("DEGRAD: mcomp  = %d \n", mcomp);
        //printf("DEGRAD: mpair  = %d \n", mpair);
        //printf("DEGRAD: n4     = %d \n", n4);
        //printf("DEGRAD: n5     = %d \n", n5);
        //printf("DEGRAD: n6     = %d \n", n6);
        //printf("DEGRAD: n7     = %d \n", n7);
        //printf("DEGRAD: n8     = %d \n", n8);
        //printf("DEGRAD: n9     = %d \n", n9);
        //printf("DEGRAD: n10    = %d \n", n10);

        // Store the truth information of the event in a file
        string photoelectron_file = "run_" + fixedLength(job) + "_" + dateb + "_" + timeb + "_photoelectrons.txt";
        string photoelectron_content = to_string(event) + "\t" + to_string(energy) + "\t" + to_string(angle) + "\t" + to_string(position) + "\t" + to_string(nclus) + "\t" + to_string(degrad_seed);
        write_text_to_position_file(photoelectron_file.c_str(), photoelectron_content.c_str());

        // For the gas gain the gain is drawn from a polya distribution
        TF1 polya = TF1("polya","([0] / [1]) *(((([1]*[1])/([2]*[2]))^(([1]*[1])/([2]*[2]))) /(TMath::Gamma((([1]*[1])/([2]*[2]))))) * ((x /[1])^((([1]*[1])/([2]*[2]))-1)) * exp(-(([1]*[1])/([2]*[2])) *(x / [1]))", 1000, 50000);
        polya.SetParameter(0, amplification_scaling); //scaling
        polya.SetParameter(1, amplification_gain); //gain
        polya.SetParameter(2, amplification_width); //width

        // Iterate over all secondary electrons of the photoelectron track and drift them to the readout
        for (Int_t iclus = 0; iclus < nclus; iclus++){
            cout << "\rGARFIELD: Electron: " << iclus + 1 << " of " << nclus << flush;
            // read in the data
            in >> x >> y >> z >> t >> n1 >> n2 >> n3;
            nlines++;
            
            // Rotate secondary electron start coordinates by defined angle
            double rot_x = x * TMath::Cos(angle) - y * TMath::Sin(angle);
            double rot_y = x * TMath::Sin(angle) + y * TMath::Cos(angle);

            if(simulation_approach == 0){
                std::unique_ptr<AvalancheMC> aval = std::make_unique<AvalancheMC>();
                aval->SetSensor(&sensor);
                aval->AvalancheElectron(rot_x / 10000., rot_y / 10000., position + (z / 10000.), t / 1000.0);
                aval->GetNumberOfElectronEndpoints();
                aval->GetElectronEndpoint(0, x1, y1, z1, t1, x2, y2, z2, t2, status);
            }
            else if(simulation_approach == 1){
                std::unique_ptr<AvalancheMicroscopic> aval = std::make_unique<AvalancheMicroscopic>();
                aval->SetSensor(&sensor);
                aval->DriftElectron(rot_x / 10000., rot_y / 10000., position + (z / 10000.), t / 1000.0, 0., 0., 0., 0.);
                aval->GetNumberOfElectronEndpoints();
                aval->GetElectronEndpoint(0, x1, y1, z1, t1, e1, x2, y2, z2, t2, e2, status);
            }
            else if(simulation_approach == 2){
                // Calculate diffusion sigma based on electron z position and diffusion coefficient
                double sigma = 10000 * el_Tdiff * TMath::Sqrt(position + (z / 10000.));
                double a = degrad_random.Gaus(0, sigma);
                double b = degrad_random.Gaus(0, sigma);

                // Get end point parameters for the electron
                x2 = (rot_x + a)/ 10000.;
                y2 = (rot_y + b)/ 10000.;
                t2 = (position + (z / 10000.))/el_vel; // Preliminary - will also depend on longitudinal diffusion
            }
            else{
                return -1;
            }

            // Draw a gas gain from the polya distribution
            double amp = polya.GetRandom();
            //cout << "\t z_start: " << position + (z / 10000.) << "\t z_stop: " << z2 << "\t electrons: " << amp << "\t status: " << status << endl; // Debug output

            // Get the pixelcoordinates oh the electron
            int posx = floor((x2 + 0.7) / pixelsize);
            int posy = floor((y2 + 0.7) / pixelsize);
            int toa = floor(t2 / 25); // 40 MHz clock
            int ftoa = 16 - floor((t2 - toa * 25) / 1.5625); // 640 MHz clock
            if(posx < 0 || posx > 255 || posy < 0 || posy > 255){
                continue;
            }
            if (hits[posx][posy] == 0){ //Counts activated pixels
                number++;
            }
            // Store the hit and the gain in a matrix
            hits[posx][posy] += (int)amp;
        }

        // Close the degrad output file and move it in the runfolder with a new name based on the eventnumber
        in.close();
        if(degrad_output == 1){
            string cmd_degrad_out = "cp DEGRAD.OUT " + degrad_dir + "/run_" + fixedLength(job) + "_data_" + fixedLength(event) + "_" + date() + "_" + time() + ".OUT";
            string err = "";
            err = exec(cmd_degrad_out.c_str());
        }

        cout << endl;

        // Write the TOS data file with the zerosupressed x, y and gain data
        f << "FEC 0\n";
        f << "Board 0\n";
        f << "Chip 1 ,Hits: " << number << "\n";
        for (int l = pixel - 1; l >= 0; l--){
            for (int m = pixel - 1; m >= 0; m--){
                if (hits[m][l] != 0){
                    f << m << " " << l << " " << hits[m][l] << "\n";
                }
            }
        }
        f.close();

        event++;
        photoelectrons++;

        // Stop if the desired amount of photoelectrons was simulated
        if(photoelectrons == nEvents){
            break;
        }
    }
}
