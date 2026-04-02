/*******************************************************************************
* Copyright (c) 2012-2013, The Microsystems Design Labratory (MDL)
* Department of Computer Science and Engineering, The Pennsylvania State University
* Exascale Computing Lab, Hewlett-Packard Company
* All rights reserved.
* 
* This source code is part of NVSim - An area, timing and power model for both 
* volatile (e.g., SRAM, DRAM) and non-volatile memory (e.g., PCRAM, STT-RAM, ReRAM, 
* SLC NAND Flash). The source code is free and you can redistribute and/or modify it
* by providing that the following conditions are met:
* 
*  1) Redistributions of source code must retain the above copyright notice,
*     this list of conditions and the following disclaimer.
* 
*  2) Redistributions in binary form must reproduce the above copyright notice,
*     this list of conditions and the following disclaimer in the documentation
*     and/or other materials provided with the distribution.
* 
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
* ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
* FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
* DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
* SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
* CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
* OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
* OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
* 
* Author list: 
*   Cong Xu	    ( Email: czx102 at psu dot edu 
*                     Website: http://www.cse.psu.edu/~czx102/ )
*   Xiangyu Dong    ( Email: xydong at cse dot psu dot edu
*                     Website: http://www.cse.psu.edu/~xydong/ )
*******************************************************************************/


#include "MemCell.h"
#include "formula.h"
#include "global.h"
#include "macros.h"
#include <math.h>
#include <iomanip>

MemCell::MemCell() {
	// TODO Auto-generated constructor stub
	memCellType         = PCRAM;
	area                = 0;
	aspectRatio         = 0;
	resistanceOn        = 0;
	resistanceOff       = 0;
	readMode            = true;
	readVoltage         = 0;
	readCurrent         = 0;
	readPower           = 0;
        wordlineBoostRatio  = 1.0;
	resetMode           = true;
	resetVoltage        = 0;
	resetCurrent        = 0;
	minSenseVoltage     = 0.08;
	resetPulse          = 0;
	resetEnergy         = 0;
	setMode             = true;
	setVoltage          = 0;
	setCurrent          = 0;
	setPulse            = 0;
	accessType          = CMOS_access;
	processNode         = 0;
	setEnergy           = 0;

	/* Optional */
	stitching         = 0;
	gateOxThicknessFactor = 2;
	widthSOIDevice = 0;
	widthAccessCMOS   = 0;
	voltageDropAccessDevice = 0;
	leakageCurrentAccessDevice = 0;
	capDRAMCell		  = 0;
	widthSRAMCellNMOS = 2.08;	/* Default NMOS width in SRAM cells is 2.08 (from CACTI) */
	widthSRAMCellPMOS = 1.23;	/* Default PMOS width in SRAM cells is 1.23 (from CACTI) */

	/*For memristors */
	readFloating = false;
	resistanceOnAtSetVoltage = 0;
	resistanceOffAtSetVoltage = 0;
	resistanceOnAtResetVoltage = 0;
	resistanceOffAtResetVoltage = 0;
	resistanceOnAtReadVoltage = 0;
	resistanceOffAtReadVoltage = 0;
	resistanceOnAtHalfReadVoltage = 0;
	resistanceOffAtHalfReadVoltage = 0;

	/* For FeDiode */
	capacitanceFeDiode        = 0;
	capacitanceFeDiodeReverse = 0;
	polarizationRemnant       = 0;
	polarizationSpontaneous   = 0;
	coerciveField             = 0;
	ferroelectricThickness    = 0;
	interlayerThickness       = 0;
	interlayerPermittivity    = 0;
	ferroelectricPermittivity = 0;
	eta                       = 0;
	ferroelectricMaterial     = "";
	electronAffinityFerroelectric = 0;
	electronAffinityInterlayer    = 0;
	effectiveMassFerroelectric    = 0;
	effectiveMassInterlayer       = 0;
	trapDepth                     = 0;
	workFunctionAnode             = 0;
	workFunctionCathode           = 0;
}

MemCell::~MemCell() {
	// TODO Auto-generated destructor stub
}

void MemCell::ReadCellFromFile(const string & inputFile)
{
	FILE *fp = fopen(inputFile.c_str(), "r");
	char line[5000];
	char tmp[5000];

	if (!fp) {
		cout << inputFile << " cannot be found!\n";
		exit(-1);
	}

	while (fscanf(fp, "%[^\n]\n", line) != EOF) {
		if (!strncmp("-MemCellType", line, strlen("-MemCellType"))) {
			sscanf(line, "-MemCellType: %s", tmp);
			if (!strcmp(tmp, "SRAM"))
				memCellType = SRAM;
			else if (!strcmp(tmp, "DRAM"))
				memCellType = DRAM;
			else if (!strcmp(tmp, "eDRAM"))
				memCellType = eDRAM;
			else if (!strcmp(tmp, "MRAM"))
				memCellType = MRAM;
			else if (!strcmp(tmp, "PCRAM"))
				memCellType = PCRAM;
			else if (!strcmp(tmp, "FBRAM"))
				memCellType = FBRAM;
			else if (!strcmp(tmp, "memristor"))
				memCellType = memristor;
			else if (!strcmp(tmp, "SLCNAND"))
				memCellType = SLCNAND;
			else if (!strcmp(tmp, "FeDiode")) {
				memCellType = FeDiode;
				voltageDropAccessDevice = 0;  /* no separate access device */
			} else
				memCellType = MLCNAND;
			continue;
		}
		if (!strncmp("-ProcessNode", line, strlen("-ProcessNode"))) {
			sscanf(line, "-ProcessNode: %d", &processNode);
			continue;
		}
		if (!strncmp("-CellArea", line, strlen("-CellArea"))) {
			sscanf(line, "-CellArea (F^2): %lf", &area);
			continue;
		}
		if (!strncmp("-CellAspectRatio", line, strlen("-CellAspectRatio"))) {
			sscanf(line, "-CellAspectRatio: %lf", &aspectRatio);
			heightInFeatureSize = sqrt(area * aspectRatio);
			widthInFeatureSize = sqrt(area / aspectRatio);
			continue;
		}

		if (!strncmp("-ResistanceOnAtSetVoltage", line, strlen("-ResistanceOnAtSetVoltage"))) {
			sscanf(line, "-ResistanceOnAtSetVoltage (ohm): %lf", &resistanceOnAtSetVoltage);
			continue;
		}
		if (!strncmp("-ResistanceOffAtSetVoltage", line, strlen("-ResistanceOffAtSetVoltage"))) {
			sscanf(line, "-ResistanceOffAtSetVoltage (ohm): %lf", &resistanceOffAtSetVoltage);
			continue;
		}
		if (!strncmp("-ResistanceOnAtResetVoltage", line, strlen("-ResistanceOnAtResetVoltage"))) {
			sscanf(line, "-ResistanceOnAtResetVoltage (ohm): %lf", &resistanceOnAtResetVoltage);
			continue;
		}
		if (!strncmp("-ResistanceOffAtResetVoltage", line, strlen("-ResistanceOffAtResetVoltage"))) {
			sscanf(line, "-ResistanceOffAtResetVoltage (ohm): %lf", &resistanceOffAtResetVoltage);
			continue;
		}
		if (!strncmp("-ResistanceOnAtReadVoltage", line, strlen("-ResistanceOnAtReadVoltage"))) {
			sscanf(line, "-ResistanceOnAtReadVoltage (ohm): %lf", &resistanceOnAtReadVoltage);
			resistanceOn = resistanceOnAtReadVoltage;
			continue;
		}
		if (!strncmp("-ResistanceOffAtReadVoltage", line, strlen("-ResistanceOffAtReadVoltage"))) {
			sscanf(line, "-ResistanceOffAtReadVoltage (ohm): %lf", &resistanceOffAtReadVoltage);
			resistanceOff = resistanceOffAtReadVoltage;
			continue;
		}
		if (!strncmp("-ResistanceOnAtHalfReadVoltage", line, strlen("-ResistanceOnAtHalfReadVoltage"))) {
			sscanf(line, "-ResistanceOnAtHalfReadVoltage (ohm): %lf", &resistanceOnAtHalfReadVoltage);
			continue;
		}
		if (!strncmp("-ResistanceOffAtHalfReadVoltage", line, strlen("-ResistanceOffAtHalfReadVoltage"))) {
			sscanf(line, "-ResistanceOffAtHalfReadVoltage (ohm): %lf", &resistanceOffAtHalfReadVoltage);
			continue;
		}
		if (!strncmp("-ResistanceOnAtHalfResetVoltage", line, strlen("-ResistanceOnAtHalfResetVoltage"))) {
			sscanf(line, "-ResistanceOnAtHalfResetVoltage (ohm): %lf", &resistanceOnAtHalfResetVoltage);
			continue;
		}

		if (!strncmp("-ResistanceOn", line, strlen("-ResistanceOn"))) {
			sscanf(line, "-ResistanceOn (ohm): %lf", &resistanceOn);
			continue;
		}
		if (!strncmp("-ResistanceOff", line, strlen("-ResistanceOff"))) {
			sscanf(line, "-ResistanceOff (ohm): %lf", &resistanceOff);
			continue;
		}
		if (!strncmp("-CapacitanceOn", line, strlen("-CapacitanceOn"))) {
			sscanf(line, "-CapacitanceOn (F): %lf", &capacitanceOn);
			continue;
		}
		if (!strncmp("-CapacitanceOff", line, strlen("-CapacitanceOff"))) {
			sscanf(line, "-CapacitanceOff (F): %lf", &capacitanceOff);
			continue;
		}

		if (!strncmp("-GateOxThicknessFactor", line, strlen("-GateOxThicknessFactor"))) {
			sscanf(line, "-GateOxThicknessFactor: %lf", &gateOxThicknessFactor);
			continue;
		}

		if (!strncmp("-SOIDeviceWidth (F)", line, strlen("-SOIDeviceWidth (F)"))) {
			sscanf(line, "-SOIDeviceWidth (F): %lf", &widthSOIDevice);
			continue;
		}

		if (!strncmp("-ReadMode", line, strlen("-ReadMode"))) {
			sscanf(line, "-ReadMode: %s", tmp);
			if (!strcmp(tmp, "voltage"))
				readMode = true;
			else
				readMode = false;
			continue;
		}
		if (!strncmp("-ReadVoltage", line, strlen("-ReadVoltage"))) {
			sscanf(line, "-ReadVoltage (V): %lf", &readVoltage);
			continue;
		}
		if (!strncmp("-ReadCurrent", line, strlen("-ReadCurrent"))) {
			sscanf(line, "-ReadCurrent (uA): %lf", &readCurrent);
			readCurrent /= 1e6;
			continue;
		}
		if (!strncmp("-ReadPower", line, strlen("-ReadPower"))) {
			sscanf(line, "-ReadPower (uW): %lf", &readPower);
			readPower /= 1e6;
			continue;
		}
		if (!strncmp("-WordlineBoostRatio", line, strlen("-WordlineBoostRatio"))) {
			sscanf(line, "-WordlineBoostRatio: %lf", &wordlineBoostRatio);
			continue;
		}
		if (!strncmp("-MinSenseVoltage", line, strlen("-MinSenseVoltage"))) {
			sscanf(line, "-MinSenseVoltage (mV): %lf", &minSenseVoltage);
			minSenseVoltage /= 1e3;
			continue;
		}


		if (!strncmp("-ResetMode", line, strlen("-ResetMode"))) {
			sscanf(line, "-ResetMode: %s", tmp);
			if (!strcmp(tmp, "voltage"))
				resetMode = true;
			else
				resetMode = false;
			continue;
		}
		if (!strncmp("-ResetVoltage", line, strlen("-ResetVoltage"))) {
			sscanf(line, "-ResetVoltage (V): %lf", &resetVoltage);
			continue;
		}
		if (!strncmp("-ResetCurrent", line, strlen("-ResetCurrent"))) {
			sscanf(line, "-ResetCurrent (uA): %lf", &resetCurrent);
			resetCurrent /= 1e6;
			continue;
		}
		if (!strncmp("-ResetVoltage", line, strlen("-ResetVoltage"))) {
			sscanf(line, "-ResetVoltage (V): %lf", &resetVoltage);
			continue;
		}
		if (!strncmp("-ResetPulse", line, strlen("-ResetPulse"))) {
			sscanf(line, "-ResetPulse (ns): %lf", &resetPulse);
			resetPulse /= 1e9;
			continue;
		}
		if (!strncmp("-ResetEnergy", line, strlen("-ResetEnergy"))) {
			sscanf(line, "-ResetEnergy (pJ): %lf", &resetEnergy);
			resetEnergy /= 1e12;
			continue;
		}

		if (!strncmp("-SetMode", line, strlen("-SetMode"))) {
			sscanf(line, "-SetMode: %s", tmp);
			if (!strcmp(tmp, "voltage"))
				setMode = true;
			else
				setMode = false;
			continue;
		}
		if (!strncmp("-SetVoltage", line, strlen("-SetVoltage"))) {
			sscanf(line, "-SetVoltage (V): %lf", &setVoltage);
			continue;
		}
		if (!strncmp("-SetCurrent", line, strlen("-SetCurrent"))) {
			sscanf(line, "-SetCurrent (uA): %lf", &setCurrent);
			setCurrent /= 1e6;
			continue;
		}
		if (!strncmp("-SetVoltage", line, strlen("-SetVoltage"))) {
			sscanf(line, "-SetVoltage (V): %lf", &setVoltage);
			continue;
		}
		if (!strncmp("-SetPulse", line, strlen("-SetPulse"))) {
			sscanf(line, "-SetPulse (ns): %lf", &setPulse);
			setPulse /= 1e9;
			continue;
		}
		if (!strncmp("-SetEnergy", line, strlen("-SetEnergy"))) {
			sscanf(line, "-SetEnergy (pJ): %lf", &setEnergy);
			setEnergy /= 1e12;
			continue;
		}

		if (!strncmp("-AccessType", line, strlen("-AccessType"))) {
			sscanf(line, "-AccessType: %s", tmp);
			if (!strcmp(tmp, "CMOS"))
				accessType = CMOS_access;
			else if (!strcmp(tmp, "BJT"))
				accessType = BJT_access;
			else if (!strcmp(tmp, "diode"))
				accessType = diode_access;
			else
				accessType = none_access;
			continue;
		}

		if (!strncmp("-AccessCMOSWidth", line, strlen("-AccessCMOSWidth"))) {
			if (accessType != CMOS_access)
				cout << "Warning: The input of CMOS access transistor width is ignored because the cell is not CMOS-accessed." << endl;
			else
				sscanf(line, "-AccessCMOSWidth (F): %lf", &widthAccessCMOS);
			continue;
		}

		if (!strncmp("-VoltageDropAccessDevice", line, strlen("-VoltageDropAccessDevice"))) {
			sscanf(line, "-VoltageDropAccessDevice (V): %lf", &voltageDropAccessDevice);
			continue;
		}

		if (!strncmp("-LeakageCurrentAccessDevice", line, strlen("-LeakageCurrentAccessDevice"))) {
			sscanf(line, "-LeakageCurrentAccessDevice (uA): %lf", &leakageCurrentAccessDevice);
			leakageCurrentAccessDevice /= 1e6;
			continue;
		}

		if (!strncmp("-DRAMCellCapacitance", line, strlen("-DRAMCellCapacitance"))) {
			if (memCellType != DRAM && memCellType != eDRAM)
				cout << "Warning: The input of DRAM cell capacitance is ignored because the memory cell is not DRAM." << endl;
			else
				sscanf(line, "-DRAMCellCapacitance (F): %lf", &capDRAMCell);
			continue;
		}

		if (!strncmp("-SRAMCellNMOSWidth", line, strlen("-SRAMCellNMOSWidth"))) {
			if (memCellType != SRAM)
				cout << "Warning: The input of SRAM cell NMOS width is ignored because the memory cell is not SRAM." << endl;
			else
				sscanf(line, "-SRAMCellNMOSWidth (F): %lf", &widthSRAMCellNMOS);
			continue;
		}

		if (!strncmp("-SRAMCellPMOSWidth", line, strlen("-SRAMCellPMOSWidth"))) {
			if (memCellType != SRAM)
				cout << "Warning: The input of SRAM cell PMOS width is ignored because the memory cell is not SRAM." << endl;
			else
				sscanf(line, "-SRAMCellPMOSWidth (F): %lf", &widthSRAMCellPMOS);
			continue;
		}


		if (!strncmp("-ReadFloating", line, strlen("-ReadFloating"))) {
			sscanf(line, "-ReadFloating: %s", tmp);
			if (!strcmp(tmp, "true"))
				readFloating = true;
			else
				readFloating = false;
			continue;
		}

		if (!strncmp("-FlashEraseVoltage (V)", line, strlen("-FlashEraseVoltage (V)"))) {
			if (memCellType != SLCNAND && memCellType != MLCNAND)
				cout << "Warning: The input of programming/erase voltage is ignored because the memory cell is not flash." << endl;
			else
				sscanf(line, "-FlashEraseVoltage (V): %lf", &flashEraseVoltage);
			continue;
		}

		if (!strncmp("-FlashProgramVoltage (V)", line, strlen("-FlashProgramVoltage (V)"))) {
			if (memCellType != SLCNAND && memCellType != MLCNAND)
				cout << "Warning: The input of programming/program voltage is ignored because the memory cell is not flash." << endl;
			else
				sscanf(line, "-FlashProgramVoltage (V): %lf", &flashProgramVoltage);
			continue;
		}

		if (!strncmp("-FlashPassVoltage (V)", line, strlen("-FlashPassVoltage (V)"))) {
			if (memCellType != SLCNAND && memCellType != MLCNAND)
				cout << "Warning: The input of pass voltage is ignored because the memory cell is not flash." << endl;
			else
				sscanf(line, "-FlashPassVoltage (V): %lf", &flashPassVoltage);
			continue;
		}

		if (!strncmp("-FlashEraseTime", line, strlen("-FlashEraseTime"))) {
			if (memCellType != SLCNAND && memCellType != MLCNAND)
				cout << "Warning: The input of erase time is ignored because the memory cell is not flash." << endl;
			else {
				sscanf(line, "-FlashEraseTime (ms): %lf", &flashEraseTime);
				flashEraseTime /= 1e3;
			}
			continue;
		}

		if (!strncmp("-FlashProgramTime", line, strlen("-FlashProgramTime"))) {
			if (memCellType != SLCNAND && memCellType != MLCNAND)
				cout << "Warning: The input of erase time is ignored because the memory cell is not flash." << endl;
			else {
				sscanf(line, "-FlashProgramTime (us): %lf", &flashProgramTime);
				flashProgramTime /= 1e6;
			}
			continue;
		}

		if (!strncmp("-GateCouplingRatio", line, strlen("-GateCouplingRatio"))) {
			if (memCellType != SLCNAND && memCellType != MLCNAND)
				cout << "Warning: The input of gate coupling ratio (GCR) is ignored because the memory cell is not flash." << endl;
			else {
				sscanf(line, "-GateCouplingRatio: %lf", &gateCouplingRatio);
			}
			continue;
		}

		/* FeDiode-specific parameters */
		if (!strncmp("-CapacitanceFeDiode", line, strlen("-CapacitanceFeDiode"))) {
			if (memCellType != FeDiode)
				cout << "Warning: -CapacitanceFeDiode is ignored because the memory cell is not FeDiode." << endl;
			else
				sscanf(line, "-CapacitanceFeDiode (F): %lf", &capacitanceFeDiode);
			continue;
		}

		if (!strncmp("-CapacitanceFeDiodeReverse", line, strlen("-CapacitanceFeDiodeReverse"))) {
			if (memCellType != FeDiode)
				cout << "Warning: -CapacitanceFeDiodeReverse is ignored because the memory cell is not FeDiode." << endl;
			else
				sscanf(line, "-CapacitanceFeDiodeReverse (F): %lf", &capacitanceFeDiodeReverse);
			continue;
		}

		if (!strncmp("-PolarizationRemnant", line, strlen("-PolarizationRemnant"))) {
			if (memCellType != FeDiode)
				cout << "Warning: -PolarizationRemnant is ignored because the memory cell is not FeDiode." << endl;
			else
				sscanf(line, "-PolarizationRemnant (uC/cm^2): %lf", &polarizationRemnant);
			continue;
		}

		if (!strncmp("-PolarizationSpontaneous", line, strlen("-PolarizationSpontaneous"))) {
			if (memCellType != FeDiode)
				cout << "Warning: -PolarizationSpontaneous is ignored because the memory cell is not FeDiode." << endl;
			else
				sscanf(line, "-PolarizationSpontaneous (uC/cm^2): %lf", &polarizationSpontaneous);
			continue;
		}

		if (!strncmp("-CoerciveField", line, strlen("-CoerciveField"))) {
			if (memCellType != FeDiode)
				cout << "Warning: -CoerciveField is ignored because the memory cell is not FeDiode." << endl;
			else
				sscanf(line, "-CoerciveField (MV/cm): %lf", &coerciveField);
			continue;
		}

		if (!strncmp("-FerroelectricThickness", line, strlen("-FerroelectricThickness"))) {
			if (memCellType != FeDiode)
				cout << "Warning: -FerroelectricThickness is ignored because the memory cell is not FeDiode." << endl;
			else
				sscanf(line, "-FerroelectricThickness (nm): %lf", &ferroelectricThickness);
			continue;
		}

		if (!strncmp("-InterlayerThickness", line, strlen("-InterlayerThickness"))) {
			if (memCellType != FeDiode)
				cout << "Warning: -InterlayerThickness is ignored because the memory cell is not FeDiode." << endl;
			else
				sscanf(line, "-InterlayerThickness (nm): %lf", &interlayerThickness);
			continue;
		}

		if (!strncmp("-InterlayerPermittivity", line, strlen("-InterlayerPermittivity"))) {
			if (memCellType != FeDiode)
				cout << "Warning: -InterlayerPermittivity is ignored because the memory cell is not FeDiode." << endl;
			else
				sscanf(line, "-InterlayerPermittivity: %lf", &interlayerPermittivity);
			continue;
		}

		if (!strncmp("-FerroelectricPermittivity", line, strlen("-FerroelectricPermittivity"))) {
			if (memCellType != FeDiode)
				cout << "Warning: -FerroelectricPermittivity is ignored because the memory cell is not FeDiode." << endl;
			else
				sscanf(line, "-FerroelectricPermittivity: %lf", &ferroelectricPermittivity);
			continue;
		}

		if (!strncmp("-Eta", line, strlen("-Eta"))) {
			if (memCellType != FeDiode)
				cout << "Warning: -Eta is ignored because the memory cell is not FeDiode." << endl;
			else
				sscanf(line, "-Eta: %lf", &eta);
			continue;
		}

		if (!strncmp("-FerroelectricMaterial", line, strlen("-FerroelectricMaterial"))) {
			if (memCellType != FeDiode)
				cout << "Warning: -FerroelectricMaterial is ignored because the memory cell is not FeDiode." << endl;
			else {
				sscanf(line, "-FerroelectricMaterial: %s", tmp);
				ferroelectricMaterial = tmp;
			}
			continue;
		}

		if (!strncmp("-ElectronAffinityFerroelectric", line, strlen("-ElectronAffinityFerroelectric"))) {
			if (memCellType != FeDiode)
				cout << "Warning: -ElectronAffinityFerroelectric is ignored because the memory cell is not FeDiode." << endl;
			else
				sscanf(line, "-ElectronAffinityFerroelectric (eV): %lf", &electronAffinityFerroelectric);
			continue;
		}

		if (!strncmp("-ElectronAffinityInterlayer", line, strlen("-ElectronAffinityInterlayer"))) {
			if (memCellType != FeDiode)
				cout << "Warning: -ElectronAffinityInterlayer is ignored because the memory cell is not FeDiode." << endl;
			else
				sscanf(line, "-ElectronAffinityInterlayer (eV): %lf", &electronAffinityInterlayer);
			continue;
		}

		if (!strncmp("-EffectiveMassFerroelectric", line, strlen("-EffectiveMassFerroelectric"))) {
			if (memCellType != FeDiode)
				cout << "Warning: -EffectiveMassFerroelectric is ignored because the memory cell is not FeDiode." << endl;
			else
				sscanf(line, "-EffectiveMassFerroelectric (m_e): %lf", &effectiveMassFerroelectric);
			continue;
		}

		if (!strncmp("-EffectiveMassInterlayer", line, strlen("-EffectiveMassInterlayer"))) {
			if (memCellType != FeDiode)
				cout << "Warning: -EffectiveMassInterlayer is ignored because the memory cell is not FeDiode." << endl;
			else
				sscanf(line, "-EffectiveMassInterlayer (m_e): %lf", &effectiveMassInterlayer);
			continue;
		}

		if (!strncmp("-TrapDepth", line, strlen("-TrapDepth"))) {
			if (memCellType != FeDiode)
				cout << "Warning: -TrapDepth is ignored because the memory cell is not FeDiode." << endl;
			else
				sscanf(line, "-TrapDepth (eV): %lf", &trapDepth);
			continue;
		}

		if (!strncmp("-WorkFunctionAnode", line, strlen("-WorkFunctionAnode"))) {
			if (memCellType != FeDiode)
				cout << "Warning: -WorkFunctionAnode is ignored because the memory cell is not FeDiode." << endl;
			else
				sscanf(line, "-WorkFunctionAnode (eV): %lf", &workFunctionAnode);
			continue;
		}

		if (!strncmp("-WorkFunctionCathode", line, strlen("-WorkFunctionCathode"))) {
			if (memCellType != FeDiode)
				cout << "Warning: -WorkFunctionCathode is ignored because the memory cell is not FeDiode." << endl;
			else
				sscanf(line, "-WorkFunctionCathode (eV): %lf", &workFunctionCathode);
			continue;
		}
	}

	/* Auto-calculate eta from Pr/Ps if not explicitly provided */
	if (memCellType == FeDiode && eta == 0 && polarizationSpontaneous > 0)
		eta = polarizationRemnant / polarizationSpontaneous;

	fclose(fp);
}


void MemCell::CellScaling(int _targetProcessNode) {
	if ((processNode > 0) && (processNode != _targetProcessNode)) {
		double scalingFactor = (double)processNode / _targetProcessNode;
		if (memCellType == PCRAM) {
			resistanceOn *= scalingFactor;
			resistanceOff *= scalingFactor;
			if (!setMode) {
				setCurrent /= scalingFactor;
			} else {
				setVoltage *= 1;
			}
			if (!resetMode) {
				resetCurrent /= scalingFactor;
			} else {
				resetVoltage *= 1;
			}
			if (accessType == diode_access) {
				capacitanceOn /= scalingFactor; //TO-DO
				capacitanceOff /= scalingFactor; //TO-DO
			}
		} else if (memCellType == MRAM){ //TO-DO: MRAM
			resistanceOn *= scalingFactor * scalingFactor;
			resistanceOff *= scalingFactor * scalingFactor;
			if (!setMode) {
				setCurrent /= scalingFactor;
			} else {
				setVoltage *= scalingFactor;
			}
			if (!resetMode) {
				resetCurrent /= scalingFactor;
			} else {
				resetVoltage *= scalingFactor;
			}
			if (accessType == diode_access) {
				capacitanceOn /= scalingFactor; //TO-DO
				capacitanceOff /= scalingFactor; //TO-DO
			}
		} else if (memCellType == memristor) { //TO-DO: memristor

		} else { //TO-DO: other RAMs

		}
		processNode = _targetProcessNode;
	}
}

double MemCell::GetMemristance(double _relativeReadVoltage) { /* Get the LRS resistance of memristor at log-linera region of I-V curve */
	if (memCellType == memristor) {
		double x1, x2, x3;  // x1: read voltage, x2: half voltage, x3: applied voltage
		if (readVoltage == 0) {
			x1 = readCurrent * resistanceOnAtReadVoltage;
		} else {
			x1 = readVoltage;
		}
		x2 = readVoltage / 2;
		x3 = _relativeReadVoltage * readVoltage;
		double y1, y2 ,y3; // y1:log(read current), y2: log(leakage current at half read voltage
		y1 = log2(x1/resistanceOnAtReadVoltage);
		y2 = log2(x2/resistanceOnAtHalfReadVoltage);
		y3 = (y2 - y1) / (x2 -x1) * x3 + (x2 * y1 - x1 * y2) / (x2 - x1);  //insertion
		return x3 / pow(2, y3);
	} else {  // not memristor, can't call the function
		cout <<"Warning[MemCell] : Try to get memristance from a non-memristor memory cell" << endl;
		return -1;
	}
}

void MemCell::CalculateWriteEnergy() {
	if (resetEnergy == 0) {
		if (resetMode) {
			if (memCellType == memristor)
				if (accessType == none_access)
					resetEnergy = fabs(resetVoltage) * (fabs(resetVoltage) - voltageDropAccessDevice) / resistanceOnAtResetVoltage * resetPulse;
				else
					resetEnergy = fabs(resetVoltage) * (fabs(resetVoltage) - voltageDropAccessDevice) / resistanceOn * resetPulse;
			else if (memCellType == PCRAM)
				resetEnergy = fabs(resetVoltage) * (fabs(resetVoltage) - voltageDropAccessDevice) / resistanceOn * resetPulse;	// PCM cells shows low resistance during most time of the switching
			else if (memCellType == FBRAM)
				resetEnergy = fabs(resetVoltage) * fabs(resetCurrent) * resetPulse;
			else
				resetEnergy = fabs(resetVoltage) * (fabs(resetVoltage) - voltageDropAccessDevice) / resistanceOn * resetPulse;
		} else {
			if (resetVoltage == 0){
				resetEnergy = tech->vdd * fabs(resetCurrent) * resetPulse; /*TO-DO consider charge pump*/
			} else {
				resetEnergy = fabs(resetVoltage) * fabs(resetCurrent) * resetPulse;
			}
			/* previous model seems to be problematic
			if (memCellType == memristor)
				if (accessType == none_access)
					resetEnergy = resetCurrent * (resetCurrent * resistanceOffAtResetVoltage + voltageDropAccessDevice) * resetPulse;
				else
					resetEnergy = resetCurrent * (resetCurrent * resistanceOff + voltageDropAccessDevice) * resetPulse;
			else if (memCellType == PCRAM)
				resetEnergy = resetCurrent * (resetCurrent * resistanceOn + voltageDropAccessDevice) * resetPulse;		// PCM cells shows low resistance during most time of the switching
			else if (memCellType == FBRAM)
				resetEnergy = fabs(resetVoltage) * fabs(resetCurrent) * resetPulse;
			else
				resetEnergy = resetCurrent * (resetCurrent * resistanceOff + voltageDropAccessDevice) * resetPulse;
		    */
		}
	}
	if (setEnergy == 0) {
		if (setMode) {
			if (memCellType == memristor)
				if (accessType == none_access)
					setEnergy = fabs(setVoltage) * (fabs(setVoltage) - voltageDropAccessDevice) / resistanceOnAtSetVoltage * setPulse;
				else
					setEnergy = fabs(setVoltage) * (fabs(setVoltage) - voltageDropAccessDevice) / resistanceOn * setPulse;
			else if (memCellType == PCRAM)
				setEnergy = fabs(setVoltage) * (fabs(setVoltage) - voltageDropAccessDevice) / resistanceOn * setPulse;			// PCM cells shows low resistance during most time of the switching
			else if (memCellType == FBRAM)
				setEnergy = fabs(setVoltage) * fabs(setCurrent) * setPulse;
			else
				setEnergy = fabs(setVoltage) * (fabs(setVoltage) - voltageDropAccessDevice) / resistanceOn * setPulse;
		} else {
			if (resetVoltage == 0){
				setEnergy = tech->vdd * fabs(setCurrent) * setPulse; /*TO-DO consider charge pump*/
			} else {
				setEnergy = fabs(setVoltage) * fabs(setCurrent) * setPulse;
			}
			/* previous model seems to be problematic
			if (memCellType == memristor)
				if (accessType == none_access)
					setEnergy = setCurrent * (setCurrent * resistanceOffAtSetVoltage + voltageDropAccessDevice) * setPulse;
				else
					setEnergy = setCurrent * (setCurrent * resistanceOff + voltageDropAccessDevice) * setPulse;
			else if (memCellType == PCRAM)
				setEnergy = setCurrent * (setCurrent * resistanceOn + voltageDropAccessDevice) * setPulse;		// PCM cells shows low resistance during most time of the switching
			else if (memCellType == FBRAM)
				setEnergy = fabs(setVoltage) * fabs(setCurrent) * setPulse;
			else
				setEnergy = setCurrent * (setCurrent * resistanceOff + voltageDropAccessDevice) * setPulse;
			*/
		}
	}
}

double MemCell::CalculateReadPower() { /* TO-DO consider charge pumped read voltage */
	if (readPower == 0) {
		if (cell->readMode) {	/* voltage-sensing */
			if (readVoltage == 0) { /* Current-in voltage sensing */
				return tech->vdd * readCurrent;
			}
			if (readCurrent == 0) { /*Voltage-divider sensing */
				double resInSerialForSenseAmp, maxBitlineCurrent;
				resInSerialForSenseAmp = sqrt(resistanceOn * resistanceOff);
				maxBitlineCurrent = (readVoltage - voltageDropAccessDevice) / (resistanceOn + resInSerialForSenseAmp);
				return tech->vdd * maxBitlineCurrent;
			}
		} else { /* current-sensing */
			double maxBitlineCurrent = (readVoltage - voltageDropAccessDevice) / resistanceOn;
			return tech->vdd * maxBitlineCurrent;
		}
	} else {
		return -1.0; /* should not call the function if read energy exists */
	}
	return -1.0;
}

/* ============================================================
 * CalculateMemoryWindow() — FeDiode-specific memory window analysis
 *
 * Maps the ferroelectric polarization state to the ON/OFF current ratio
 * for an MFIS two-terminal crossbar device.
 *
 * Ferroelectric physics (from Toprasertpong et al., IEEE TED 2022,
 * as adapted in Fe-NVSim):
 *   Eq. 8  — Ideal MFS window (utility function)
 *   Eq. 9  — MWminor for MFIS under voltage division
 *   Eq. 10 — alphaV, parameterizes partial switching at Vm
 *
 * Transport mapping (FeDiode-specific, NOT ΔVth FeFET model):
 *   The ferroelectric MW modulates the junction potential.
 *   Effective polarization-induced junction shift:
 *     Vpol = 2*Pr / (ε0 * (εFE/tFE + εIL/tIL))
 *   This is the MFIS capacitor voltage divider at full switching —
 *   independent of the write voltage, governed only by Pr and geometry.
 *   Upper bound on ION/IOFF = exp(q*Vpol / kT)  (ideal thermionic transport)
 *   Operational ION/IOFF = exp(q*Vpol / (n_eff*kT))  (n_eff calibrated to measured)
 *
 * Validation (20% criterion):
 *   MW_MFIS(Vwrite) must exceed the minimum MW needed for the measured
 *   ON/OFF ratio with ideal transport: MW_min = kT*ln(ION/IOFF).
 *   If not, print warning and suggest corrected Pr or Ec.
 * ============================================================ */
void MemCell::CalculateMemoryWindow()
{
	if (memCellType != FeDiode) {
		cout << "[MemCell] Warning: CalculateMemoryWindow() called on non-FeDiode cell. Skipped." << endl;
		return;
	}

	/* ---- Physical constants ---- */
	const double eps0  = 8.854e-12;  /* F/m */
	const double kT_q  = 0.02585;    /* V at 300 K */

	/* ---- Unit conversions from cell-file storage to SI ---- */
	/* polarizationRemnant:      uC/cm^2  ->  C/m^2  (×0.01) */
	/* polarizationSpontaneous:  uC/cm^2  ->  C/m^2  (×0.01) */
	/* coerciveField:            MV/cm    ->  V/m    (×1e8)  */
	/* ferroelectricThickness:   nm       ->  m      (×1e-9) */
	/* interlayerThickness:      nm       ->  m      (×1e-9) */
	double Pr  = polarizationRemnant     * 1e-2;   /* C/m^2 */
	double Ps  = polarizationSpontaneous * 1e-2;   /* C/m^2 */
	double Ec  = coerciveField           * 1e8;    /* V/m   */
	double tFE = ferroelectricThickness  * 1e-9;   /* m     */
	double tIL = interlayerThickness     * 1e-9;   /* m     */
	double eFE = ferroelectricPermittivity;        /* dimensionless */
	double eIL = interlayerPermittivity;           /* dimensionless */

	if (Pr == 0 || Ps == 0 || Ec == 0 || tFE == 0 || eFE == 0) {
		cout << "[MemCell] Warning: FeDiode memory window skipped — ferroelectric parameters incomplete." << endl;
		cout << "          Set PolarizationRemnant, PolarizationSpontaneous, CoerciveField," << endl;
		cout << "          FerroelectricThickness, and FerroelectricPermittivity in the .cell file." << endl;
		return;
	}

	/* η: squareness; auto-calculated from Pr/Ps if -Eta not given */
	double eta_val = (eta > 0) ? eta : Pr / Ps;
	if (eta_val <= 0 || eta_val > 1) {
		cout << "[MemCell] Warning: eta = " << eta_val << " out of (0,1]. Resetting to Pr/Ps = " << Pr/Ps << endl;
		eta_val = Pr / Ps;
	}

	/* ================================================================
	 * Eq. 8 — Ideal MFS memory window (utility; used as MW0 in Eq. 9)
	 *
	 * Branch 1 (small Pr):  MW = 2*Pr / (eFE*eps0)
	 * Branch 2 (large Pr):  MW = 2*Ec*tFE / (1 + eFE*eps0*Ec*tanh(η)/Pr) * (2 - 1/η)
	 *
	 * Boundary: Pr vs. eFE*eps0*Ec.  Paper says "Pr << eFE*eps0*Ec" for branch 1
	 * and "Pr > eFE*eps0*Ec*tanh(η)" for branch 2.  We use the tanh-weighted
	 * boundary as the crossover.
	 * ================================================================ */
	double eFE_eps0_Ec = eFE * eps0 * Ec;   /* C/m^2: characteristic polarization scale */
	double mw_mfs;
	if (Pr < eFE_eps0_Ec) {
		/* Small-Pr limit */
		mw_mfs = 2.0 * Pr / (eFE * eps0);
	} else {
		/* Large-Pr limit (also used when Pr is between the two branches, per paper note) */
		double factor = 1.0 + eFE * eps0 * Ec * tanh(eta_val) / Pr;
		mw_mfs = 2.0 * Ec * tFE / factor * (2.0 - 1.0 / eta_val);
	}

	/* ================================================================
	 * MFIS voltage division: Vm (voltage across FE layer) from Vg (applied)
	 *
	 *   Vm = Vg * (tFE/eFE) / (tFE/eFE + tIL/eIL)
	 *
	 * Physical picture: series capacitors C_FE = eFE*eps0/tFE and C_IL = eIL*eps0/tIL.
	 * The FE layer sees a reduced voltage relative to the total applied Vg.
	 * ================================================================ */
	double Vg  = fabs(resetVoltage);   /* write voltage magnitude */
	double rFE = tFE / eFE;            /* "electrical thickness" of FE layer */
	double rIL = (tIL > 0 && eIL > 0) ? tIL / eIL : 0.0;
	double Vm  = (rIL > 0.0) ? Vg * rFE / (rFE + rIL) : Vg;

	/* ================================================================
	 * θ- (theta-minus): shape parameter for MWminor formula (Eq. 9)
	 *   θ- ≈ (1/2) * (1 + (3/5) * eFE*eps0*Ec / (Ps*η))
	 * ================================================================ */
	double theta_m = 0.5 * (1.0 + 0.6 * eFE_eps0_Ec / (Ps * eta_val));

	/* ================================================================
	 * Transport mapping (FeDiode-specific)
	 *
	 * The remanent polarization Pr creates a bound surface charge ±2*Pr
	 * at the FE/IL interface.  For an MFIS capacitor stack at rest
	 * (no applied voltage after write), this charge induces a potential:
	 *
	 *   Vpol = 2*Pr / (eps0 * (eFE/tFE + eIL/tIL))
	 *        = 2*Pr / (eps0 / rFE + eps0 / rIL)   [series dielectrics]
	 *
	 * This Vpol shifts the diode I-V curve between the two polarization
	 * states, giving: ION/IOFF = exp(q*Vpol / (n_eff*kT))
	 *
	 * n_eff is back-calculated from the measured resistance ratio and Vpol.
	 * It encapsulates all non-ideal effects: interface traps, partial
	 * depolarization, non-thermionic transport.
	 * ================================================================ */
	double C_stack_inv = rFE / eps0 + ((rIL > 0) ? rIL / eps0 : 0.0);  /* 1/C', m^2/F */
	double Vpol = (C_stack_inv > 0) ? 2.0 * Pr * C_stack_inv : 0.0;    /* V */

	double measured_ratio = (resistanceOn > 0) ? (resistanceOff / resistanceOn) : 4000.0;
	double n_eff_vpol = 1.0;
	if (Vpol > 0 && measured_ratio > 1)
		n_eff_vpol = Vpol / (kT_q * log(measured_ratio));

	/* Minimum MW needed for measured ratio at ideal transport (n=1) */
	double mw_min_ideal = kT_q * log(measured_ratio);   /* V */

	/* ================================================================
	 * Sweep Vm from 0.1*Vm to Vm in 10 steps
	 * At each step, compute alphaV (Eq. 10) and MWminor (Eq. 9)
	 * ================================================================ */
	double mw_minor_at_vwrite = 0.0;
	int    numSteps = 10;
	double denom_tanh_Pr = tanh(theta_m * Pr / Ps);   /* denominator, constant */
	double vm_sweep[16];
	double alpha_sweep[16];
	double mw_sweep[16];

	for (int i = 1; i <= numSteps; i++) {
		double Vm_i = Vm * (double)i / numSteps;

		/* Eq. 10: alphaV — effective polarization amplitude at partial voltage Vm_i */
		double E_norm_plus  = eta_val * (Vm_i / tFE + Ec) / Ec;
		double E_norm_minus = eta_val * (Vm_i / tFE - Ec) / Ec;
		double alphaV = 0.5 * Ps * (tanh(E_norm_plus) - tanh(E_norm_minus));

		/* Eq. 9: MWminor */
		double mw_minor = 0.0;
		if (fabs(denom_tanh_Pr) > 1e-15 && alphaV > 0.0) {
			double ratio_tanh = tanh(theta_m * alphaV / Ps) / denom_tanh_Pr;
			mw_minor = mw_mfs * (1.0 - ratio_tanh)
			         - (2.0 * tFE * alphaV / (eFE * eps0)) * (1.0 - (Pr / alphaV) * ratio_tanh);
		}

		/* Clamp to zero if negative (physically: write voltage below Ec, no switching) */
		if (mw_minor < 0.0) mw_minor = 0.0;

		vm_sweep[i-1] = Vm_i;
		alpha_sweep[i-1] = alphaV;
		mw_sweep[i-1] = mw_minor;

		if (i == numSteps)
			mw_minor_at_vwrite = mw_minor;
	}

	/* Calibrate n_eff at operating voltage so ON/OFF_op tracks the measured anchor.
	 * This keeps the MW sweep shape from Eq.9/10 while using a physically interpretable
	 * transport coupling factor fitted at the operating point. */
	double n_eff = n_eff_vpol;
	if (mw_minor_at_vwrite > 0 && measured_ratio > 1)
		n_eff = mw_minor_at_vwrite / (kT_q * log(measured_ratio));

	/* ================================================================
	 * Print header
	 * ================================================================ */
	cout << endl;
	cout << "============================================================" << endl;
	cout << " FeDiode Memory Window Analysis" << endl;
	if (!ferroelectricMaterial.empty())
		cout << " Material        : " << ferroelectricMaterial << endl;
	cout << fixed << setprecision(3);
	cout << " Pr = "  << polarizationRemnant      << " uC/cm^2"
	     << "   Ps = " << polarizationSpontaneous  << " uC/cm^2"
	     << "   eta = " << eta_val << endl;
	cout << " Ec = "  << coerciveField            << " MV/cm"
	     << "   tFE = " << ferroelectricThickness  << " nm"
	     << "   eFE = " << eFE << endl;
	if (tIL > 0)
		cout << " tIL = " << interlayerThickness << " nm"
		     << "   eIL = " << eIL << endl;
	cout << " Vwrite = " << Vg << " V"
	     << "   Vm (FE layer at Vwrite) = " << Vm << " V" << endl;
	/* Print transport parameters if provided */
	if (workFunctionAnode > 0 || workFunctionCathode > 0)
		cout << " WF anode = " << workFunctionAnode << " eV"
		     << "   WF cathode = " << workFunctionCathode << " eV" << endl;
	if (electronAffinityFerroelectric > 0 || electronAffinityInterlayer > 0)
		cout << " chi_FE = " << electronAffinityFerroelectric << " eV"
		     << "   chi_IL = " << electronAffinityInterlayer << " eV" << endl;
	if (effectiveMassFerroelectric > 0 || effectiveMassInterlayer > 0)
		cout << " m*_FE = " << effectiveMassFerroelectric << " m_e"
		     << "   m*_IL = " << effectiveMassInterlayer << " m_e" << endl;
	if (trapDepth > 0)
		cout << " Trap depth = " << trapDepth << " eV" << endl;
	cout << "------------------------------------------------------------" << endl;
	cout << " MW_MFS  (ideal MFS,  Eq.8)          = " << mw_mfs << " V" << endl;
	cout << " Vpol    (polarization junction shift)= " << Vpol   << " V" << endl;
	cout << " ION/IOFF upper bound (n_eff=1)       = ";
	cout << scientific << setprecision(3) << exp(Vpol / kT_q) << endl;
	cout << fixed << setprecision(3);
	cout << " Measured R_off/R_on                  = " << measured_ratio << endl;
	cout << " n_eff (from Vpol anchor)             = " << n_eff_vpol << endl;
	cout << " n_eff (operating-point calibrated)   = " << n_eff << endl;
	cout << "------------------------------------------------------------" << endl;

	cout << endl;
	cout << "  Vm (V)   alphaV (C/m^2)  MW_minor (V)  ION/IOFF_op  ION/IOFF_max" << endl;
	cout << "  ------   -------------   ------------  -----------  ------------" << endl;

	for (int i = 0; i < numSteps; i++) {
		double ioff_op  = (n_eff > 0) ? exp(mw_sweep[i] / (n_eff * kT_q)) : 1.0;
		double ioff_max = exp(mw_sweep[i] / kT_q);

		cout << "  " << fixed << setprecision(4) << vm_sweep[i]
		     << "     " << scientific << setprecision(3) << alpha_sweep[i]
		     << "       " << fixed << setprecision(4) << mw_sweep[i]
		     << "     " << scientific << setprecision(3) << ioff_op
		     << "   " << scientific << setprecision(3) << ioff_max << endl;
	}

	/* ================================================================
	 * Validation: does MW_MFIS(Vwrite) >= MW_min needed for measured ratio?
	 * The 20% tolerance applies to the ratio MW_minor_at_vwrite / mw_min_ideal.
	 * ================================================================ */
	cout << endl;
	cout << "------------------------------------------------------------" << endl;
	cout << " Validation checkpoint" << endl;
	cout << "   MW_minor at Vwrite                = " << fixed << setprecision(4) << mw_minor_at_vwrite << " V" << endl;
	cout << "   MW_min needed (ideal transport)   = " << mw_min_ideal << " V  (for ION/IOFF = " << measured_ratio << ")" << endl;
	cout << "   Margin  MW_minor / MW_min          = " << fixed << setprecision(2) << mw_minor_at_vwrite / mw_min_ideal << "x" << endl;

	bool pass = (mw_minor_at_vwrite >= mw_min_ideal * 0.80);  /* 20% tolerance */
	if (pass) {
		cout << " [PASS] Ferroelectric parameters are consistent with measured ON/OFF." << endl;
		cout << "        n_eff = " << fixed << setprecision(2) << n_eff;
		if (n_eff <= 5.0)
			cout << " (thermionic / Schottky-like transport)" << endl;
		else if (n_eff <= 15.0)
			cout << " (trap-assisted or tunneling transport, physically plausible)" << endl;
		else
			cout << " (very low FE-junction coupling — check interface quality)" << endl;
	} else {
		cout << " [WARN] MW_minor at Vwrite is below the ideal-transport minimum." << endl;
		cout << "        The measured ON/OFF cannot be explained without supplementary" << endl;
		cout << "        conduction mechanisms, or the FE parameters need refinement." << endl;
		cout << endl;
		/* Suggest what Ec is needed to achieve MW_minor = mw_min_ideal */
		/* In the large-Pr regime: MW ≈ 2*Ec*tFE * (2 - 1/η) / (1 + eFE*eps0*Ec*tanh(η)/Pr) */
		/* For small eFE*eps0*Ec*tanh(η)/Pr: MW ≈ 2*Ec*tFE*(2-1/η) → Ec ≈ MW / (2*tFE*(2-1/η)) */
		double shape = 2.0 - 1.0 / eta_val;
		if (shape > 0 && tFE > 0) {
			double Ec_needed = mw_min_ideal / (2.0 * tFE * shape);
			cout << "        Suggested Ec to achieve MW_min with current Pr/η:" << endl;
			cout << "          Ec_needed ≈ " << fixed << setprecision(3) << Ec_needed * 1e-8 << " MV/cm" << endl;
		}
		/* Suggest Pr needed: for large-Pr regime, MW saturates; for small-Pr: MW = 2*Pr/(eFE*eps0) */
		/* → Pr_needed = mw_min_ideal * eFE * eps0 / 2 */
		double Pr_needed_SI = mw_min_ideal * eFE * eps0 / 2.0;
		cout << "        Suggested Pr (small-Pr regime limit):" << endl;
		cout << "          Pr_needed ≈ " << fixed << setprecision(2) << Pr_needed_SI * 100.0 << " uC/cm^2" << endl;
	}
	cout << "============================================================" << endl;
	cout << endl;
}

void MemCell::PrintCell()
{
	switch (memCellType) {
	case SRAM:
		cout << "Memory Cell: SRAM" << endl;
		break;
	case DRAM:
		cout << "Memory Cell: DRAM" << endl;
		break;
	case eDRAM:
		cout << "Memory Cell: Embedded DRAM" << endl;
		break;
	case MRAM:
		cout << "Memory Cell: MRAM (Magnetoresistive)" << endl;
		break;
	case PCRAM:
		cout << "Memory Cell: PCRAM (Phase-Change)" << endl;
		break;
	case memristor:
		cout << "Memory Cell: RRAM (Memristor)" << endl;
		break;
	case FBRAM:
		cout << "Memory Cell: FBRAM (Floating Body)" <<endl;
		break;
	case SLCNAND:
		cout << "Memory Cell: Single-Level Cell NAND Flash" << endl;
		break;
	case MLCNAND:
		cout << "Memory Cell: Multi-Level Cell NAND Flash" << endl;
		break;
	case FeDiode:
		cout << "Memory Cell: FeDiode (Two-terminal Ferroelectric Diode Crossbar)" << endl;
		break;
	default:
		cout << "Memory Cell: Unknown" << endl;
	}
	cout << "Cell Area (F^2)    : " << area << " (" << heightInFeatureSize << "Fx" << widthInFeatureSize << "F)" << endl;
	cout << "Cell Aspect Ratio  : " << aspectRatio << endl;

	if (memCellType == PCRAM || memCellType == MRAM || memCellType == memristor || memCellType == FBRAM) {
		if (resistanceOn < 1e3 )
			cout << "Cell Turned-On Resistance : " << resistanceOn << "ohm" << endl;
		else if (resistanceOn < 1e6)
			cout << "Cell Turned-On Resistance : " << resistanceOn / 1e3 << "Kohm" << endl;
		else
			cout << "Cell Turned-On Resistance : " << resistanceOn / 1e6 << "Mohm" << endl;
		if (resistanceOff < 1e3 )
			cout << "Cell Turned-Off Resistance: "<< resistanceOff << "ohm" << endl;
		else if (resistanceOff < 1e6)
			cout << "Cell Turned-Off Resistance: "<< resistanceOff / 1e3 << "Kohm" << endl;
		else
			cout << "Cell Turned-Off Resistance: "<< resistanceOff / 1e6 << "Mohm" << endl;

		if (readMode) {
			cout << "Read Mode: Voltage-Sensing" << endl;
			if (readCurrent > 0)
				cout << "  - Read Current: " << readCurrent * 1e6 << "uA" << endl;
			if (readVoltage > 0)
				cout << "  - Read Voltage: " << readVoltage << "V" << endl;
		} else {
			cout << "Read Mode: Current-Sensing" << endl;
			if (readCurrent > 0)
				cout << "  - Read Current: " << readCurrent * 1e6 << "uA" << endl;
			if (readVoltage > 0)
				cout << "  - Read Voltage: " << readVoltage << "V" << endl;
		}

		if (resetMode) {
			cout << "Reset Mode: Voltage" << endl;
			cout << "  - Reset Voltage: " << resetVoltage << "V" << endl;
		} else {
			cout << "Reset Mode: Current" << endl;
			cout << "  - Reset Current: " << resetCurrent * 1e6 << "uA" << endl;
		}
		cout << "  - Reset Pulse: " << TO_SECOND(resetPulse) << endl;

		if (setMode) {
			cout << "Set Mode: Voltage" << endl;
			cout << "  - Set Voltage: " << setVoltage << "V" << endl;
		} else {
			cout << "Set Mode: Current" << endl;
			cout << "  - Set Current: " << setCurrent * 1e6 << "uA" << endl;
		}
		cout << "  - Set Pulse: " << TO_SECOND(setPulse) << endl;

		switch (accessType) {
		case CMOS_access:
			cout << "Access Type: CMOS" << endl;
			break;
		case BJT_access:
			cout << "Access Type: BJT" << endl;
			break;
		case diode_access:
			cout << "Access Type: Diode" << endl;
			break;
		default:
			cout << "Access Type: None Access Device" << endl;
		}
	} else if (memCellType == SRAM) {
		cout << "SRAM Cell Access Transistor Width: " << widthAccessCMOS << "F" << endl;
		cout << "SRAM Cell NMOS Width: " << widthSRAMCellNMOS << "F" << endl;
		cout << "SRAM Cell PMOS Width: " << widthSRAMCellPMOS << "F" << endl;
	} else if (memCellType == SLCNAND) {
		cout << "Pass Voltage       : " << flashPassVoltage << "V" << endl;
		cout << "Programming Voltage: " << flashProgramVoltage << "V" << endl;
		cout << "Erase Voltage      : " << flashEraseVoltage << "V" << endl;
		cout << "Programming Time   : " << TO_SECOND(flashProgramTime) << endl;
		cout << "Erase Time         : " << TO_SECOND(flashEraseTime) << endl;
		cout << "Gate Coupling Ratio: " << gateCouplingRatio << endl;
	}
}
