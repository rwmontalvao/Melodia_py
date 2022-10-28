// Copyright 2021 KU Leuven.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
// Author: Rinaldo Wander Montalvão, PhD
//
module pdb;

/** It assumes the following PDB format ********************************************** 
 *  1 -  6        Record name     "ATOM  "
 *  7 - 11        Integer         serial        Atom serial number.
 * 13 - 16        Atom            name          Atom name.
 * 17             Character       altLoc        Alternate location indicator.
 * 18 - 20        Residue name    resName       Residue name.
 * 22             Character       chainID       Chain identifier.
 * 23 - 26        Integer         resSeq        Residue sequence number.
 * 27             AChar           iCode         Code for insertion of residues.
 * 31 - 38        Real(8.3)       x             Orthogonal coordinates for X
 * 39 - 46        Real(8.3)       y             Orthogonal coordinates for Y
 * 47 - 54        Real(8.3)       z             Orthogonal coordinates for Z
 * 55 - 60        Real(6.2)       occupancy     Occupancy.
 * 61 - 66        Real(6.2)       tempFactor    Temperature factor.
 * 73 - 76        LString(4)      segID         Segment identifier, left-justified.
 * 77 - 78        LString(2)      element       Element symbol, right-justified.
 * 79 - 80        LString(2)      charge        Charge on the atom.
 *************************************************************************************/

import geometry;

import std.algorithm;
import std.conv;
import std.math;
import std.parallelism;
import std.stdio;
import std.string;

import progress;

class Atom
{
	int serial; // Atom serial number.

	string name; // Atom name.

	char[1] altLoc; // Alternate location indicator.

	string resName; // Residue name.

	char[1] chainID; // Chain identifier.

	char[1] iCode; // Code for insertion of residues.

	double x; // Orthogonal coordinates for X.
	double y; // Orthogonal coordinates for Y.
	double z; // Orthogonal coordinates for Z.

	double occupancy; // Occupancy.
	double tempFactor; // Temperature factor.

	string segID; // Segment identifier, left-justified.
	string element; // Element symbol, right-justified.
	string charge; // Charge on the atom.

	this(in string name)
	{
		this.name = name;
	}
}

class Residue
{

	int resSeq; // Residue sequence number.

	string name; // Residue name.

	Atom[string] atoms; // Residue's atoms.

	// Geometry
	double curvature;
	double torsion;
	double writhing;
	double arc_length;

	// Torsion angles
	double phi;
	double psi;

	this(in int resSeq)
	{
		this.resSeq = resSeq;
	}
}

class Model
{

	int model_number;

	Residue[int] residues; // Molecule's residues; type of Associative arrays

	this(in int model_number)
	{
		this.model_number = model_number;
	}
}

class PDB
{
	string ID; // PDB id

	Model[int] models;

	// Create jobs to process a range of models.
	private struct Range
	{

		int a, b; // Models from a to b 

		// Job contructor
		this(in int a, in int b)
		{
			this.a = to!int(a);
			this.b = to!int(b);
		}

		void xgeo(ref Model[int] models)
		{

			for (int model_number = this.a; model_number <= this.b; model_number++)
			{

				//writefln("Calculating model %s", model_number);
				//stdout.flush();
				auto keys = sort!("a < b")(models[model_number].residues.keys);

				const int f_residue = 0;
				const int l_residue = to!int(keys.length) - 1;

				// Calculate dihedral angles
				double[4][3] list_phi;
				double[4][3] list_psi;

				foreach (i; 0 .. keys.length)
				{

					auto res_num = keys[i];

					double phi = 0.0;
					double psi = 0.0;

					Residue res = models[model_number].residues[res_num];

					list_phi[0][1] = res.atoms["N"].x;
					list_phi[1][1] = res.atoms["N"].y;
					list_phi[2][1] = res.atoms["N"].z;
					list_psi[0][0] = res.atoms["N"].x;
					list_psi[1][0] = res.atoms["N"].y;
					list_psi[2][0] = res.atoms["N"].z;

					list_phi[0][2] = res.atoms["CA"].x;
					list_phi[1][2] = res.atoms["CA"].y;
					list_phi[2][2] = res.atoms["CA"].z;
					list_psi[0][1] = res.atoms["CA"].x;
					list_psi[1][1] = res.atoms["CA"].y;
					list_psi[2][1] = res.atoms["CA"].z;

					list_phi[0][3] = res.atoms["C"].x;
					list_phi[1][3] = res.atoms["C"].y;
					list_phi[2][3] = res.atoms["C"].z;
					list_psi[0][2] = res.atoms["C"].x;
					list_psi[1][2] = res.atoms["C"].y;
					list_psi[2][2] = res.atoms["C"].z;

					if (i == f_residue)
					{
						Residue next_res = models[model_number].residues[keys[i + 1]];

						list_psi[0][3] = next_res.atoms["N"].x;
						list_psi[1][3] = next_res.atoms["N"].y;
						list_psi[2][3] = next_res.atoms["N"].z;

						psi = dihed(list_psi);
						phi = 0.0;
					}
					else if (i == l_residue)
					{
						Residue prev_res = models[model_number].residues[keys[i - 1]];

						list_phi[0][0] = prev_res.atoms["C"].x;
						list_phi[1][0] = prev_res.atoms["C"].y;
						list_phi[2][0] = prev_res.atoms["C"].z;

						psi = 0.0;
						phi = dihed(list_phi);
					}
					else
					{
						Residue next_res = models[model_number].residues[keys[i + 1]];
						Residue prev_res = models[model_number].residues[keys[i - 1]];

						list_psi[0][3] = next_res.atoms["N"].x;
						list_psi[1][3] = next_res.atoms["N"].y;
						list_psi[2][3] = next_res.atoms["N"].z;

						list_phi[0][0] = prev_res.atoms["C"].x;
						list_phi[1][0] = prev_res.atoms["C"].y;
						list_phi[2][0] = prev_res.atoms["C"].z;

						psi = dihed(list_psi);
						phi = dihed(list_phi);
					}
					models[model_number].residues[res_num].psi = psi;
					models[model_number].residues[res_num].phi = phi;

					//writefln("%s %s %s %s", res_num, models[model_number].residues[res_num].name, phi, psi);
				}

				// Differential Geometry
				immutable int n = to!int(models[model_number].residues.length);

				double[] t = new double[](n);
				double[] x = new double[](n);
				double[] y = new double[](n);
				double[] z = new double[](n);

				double[][] x2 = new double[][](n, 3);
				double[][] y2 = new double[][](n, 3);
				double[][] z2 = new double[][](n, 3);

				int i = 0, last;
				foreach (res_num; keys)
				{
					Residue res = models[model_number].residues[res_num];

					t[i] = res.resSeq;
					x[i] = res.atoms["CA"].x;
					y[i] = res.atoms["CA"].y;
					z[i] = res.atoms["CA"].z;
					last = i;
					i++;
				}

				spline(t, x, x2);
				spline(t, y, y2);
				spline(t, z, z2);

				// temp
				if (model_number == 1)
				{
					File out_file = File("model1.csv", "w");
					auto pos = t[0];
					while (pos <= t[last])
					{
						auto xs = splint(t, x, x2, pos);
						auto ys = splint(t, y, y2, pos);
						auto zs = splint(t, z, z2, pos);
						out_file.writefln("%f,%f,%f,%f", pos, xs, ys, zs);
						pos += 0.01;
					}
					out_file.close();
				}
				// end temp

				double[3] dx, dy, dz;

				i = 0;
				foreach (res_num; keys)
				{
					Residue res = models[model_number].residues[res_num];

					//writeln(res_num);

					if (i > 0 && i < last)
					{
						nderiv(t, x, x2, res.resSeq, dx);
						nderiv(t, y, y2, res.resSeq, dy);
						nderiv(t, z, z2, res.resSeq, dz);

						// Compute curvature
						double[3] v1, v2;

						v1[0] = dx[0];
						v1[1] = dy[0];
						v1[2] = dz[0];

						v2[0] = dx[1];
						v2[1] = dy[1];
						v2[2] = dz[1];

						auto rs = cross(v1, v2);
						auto r1 = dot(rs, rs);
						auto r2 = dot(v1, v1);

						res.curvature = sqrt(r1) / pow(sqrt(r2), 3);

						//Compute torsion
						immutable auto det = -dx[2] * dy[1] * dz[0]
							+ dx[1] * dy[2] * dz[0]
							+ dx[2] * dy[0] * dz[1]
							- dx[0] * dy[2] * dz[1]
							- dx[1] * dy[0] * dz[2]
							+ dx[0] * dy[1] * dz[2];

						res.torsion = det / r1;

						// Compute arc lenght
						double[3] dist;

						double sum = 0.0;

						double j = t[i - 1];
						do
						{
							dist[0] = splint(t, x, x2, j + 0.1) - splint(t, x, x2, j);
							dist[1] = splint(t, y, y2, j + 0.1) - splint(t, y, y2, j);
							dist[2] = splint(t, z, z2, j + 0.1) - splint(t, z, z2, j);

							sum += sqrt(dot(dist, dist));

							j += 0.1;
						}
						while (j < t[i + 1]);

						res.arc_length = sum;

						// Compute the writhing number
						int start, stop;

						if (i == 1)
						{
							start = i - 1;
							stop = i + 3;
						}
						else if (i == (last - 1))
						{
							start = i - 3;
							stop = i + 1;
						}
						else
						{
							start = i - 2;
							stop = i + 2;
						}

						start = to!int(t[start]);
						stop = to!int(t[stop]);

						double[3] rij, ri1j, rij1, ri1j1, rjj1, rii1;
						sum = 0.0;
						for (int ii = start; ii <= stop - 3; ii++)
						{
							for (int jj = ii + 2; jj <= stop - 1; jj++)
							{
								rij[0] = splint(t, x, x2, to!double(jj)) - splint(t, x, x2, to!double(
										ii));
								rij[1] = splint(t, y, y2, to!double(jj)) - splint(t, y, y2, to!double(
										ii));
								rij[2] = splint(t, z, z2, to!double(jj)) - splint(t, z, z2, to!double(
										ii));

								ri1j[0] = splint(t, x, x2, to!double(jj)) - splint(t, x, x2, to!double(
										ii + 1));
								ri1j[1] = splint(t, y, y2, to!double(jj)) - splint(t, y, y2, to!double(
										ii + 1));
								ri1j[2] = splint(t, z, z2, to!double(jj)) - splint(t, z, z2, to!double(
										ii + 1));

								rij1[0] = splint(t, x, x2, to!double(jj + 1)) - splint(t, x, x2, to!double(
										ii));
								rij1[1] = splint(t, y, y2, to!double(jj + 1)) - splint(t, y, y2, to!double(
										ii));
								rij1[2] = splint(t, z, z2, to!double(jj + 1)) - splint(t, z, z2, to!double(
										ii));

								ri1j1[0] = splint(t, x, x2, to!double(jj + 1)) - splint(t, x, x2, to!double(
										ii + 1));
								ri1j1[1] = splint(t, y, y2, to!double(jj + 1)) - splint(t, y, y2, to!double(
										ii + 1));
								ri1j1[2] = splint(t, z, z2, to!double(jj + 1)) - splint(t, z, z2, to!double(
										ii + 1));

								rjj1[0] = splint(t, x, x2, to!double(jj + 1)) - splint(t, x, x2, to!double(
										jj));
								rjj1[1] = splint(t, y, y2, to!double(jj + 1)) - splint(t, y, y2, to!double(
										jj));
								rjj1[2] = splint(t, z, z2, to!double(jj + 1)) - splint(t, z, z2, to!double(
										jj));

								rii1[0] = splint(t, x, x2, to!double(ii + 1)) - splint(t, x, x2, to!double(
										ii));
								rii1[1] = splint(t, y, y2, to!double(ii + 1)) - splint(t, y, y2, to!double(
										ii));
								rii1[2] = splint(t, z, z2, to!double(ii + 1)) - splint(t, z, z2, to!double(
										ii));

								auto aij = div(cross(rij, rij1), norm(cross(rij, rij1)));
								auto bij = div(cross(rij1, ri1j1), norm(cross(rij1, ri1j1)));
								auto cij = div(cross(ri1j1, ri1j), norm(cross(ri1j1, ri1j)));
								auto dij = div(cross(ri1j, rij), norm(cross(ri1j, rij)));

								immutable double Omegaij = (asin(dot(aij, bij))
										+ asin(dot(bij, cij))
										+ asin(dot(cij, dij))
										+ asin(dot(dij, aij))) * sgn(dot(cross(rjj1, rii1), rij1));

								sum += Omegaij / (4.0 * PI);
							}

						}

						res.writhing = 2.0 * sum;

					}
					else
					{
						res.curvature = 0.0;
						res.torsion = 0.0;
						res.arc_length = 0.0;
						res.writhing = 0.0;
					}
					i++;
				}
			}
		}
	}

	// Class constructor
	this(in string ID, in string file_name, in int nCPUs)
	{
		// List of known residues
		const auto residue_codes = [
			"CYS",
			"ASP",
			"ALA",
			"GLU",
			"PHE",
			"GLY",
			"HIS",
			"ILE",
			"LYS",
			"LEU",
			"MET",
			"ASN",
			"PRO",
			"GLN",
			"ARG",
			"SER",
			"THR",
			"VAL",
			"TRP",
			"TYR",
			"MSE"
		];

		const auto atom_names = ["N", "CA", "C"];

		int last_residue = -999;

		int model_number = 1;

		bool model_new = false;

		this.ID = ID;

		int count = 0;
		foreach (line; File(file_name).byLine())
			count++;

		Bar bar = new Bar();
		bar.message = { return "Loading"; };
		bar.suffix = { return bar.percent.to!string ~ "%"; };
		bar.max = count;

		//Parse PDB atoms
		foreach (line; File(file_name).byLine())
		{
			bar.next();

			// Select model
			if (line.indexOf("MODEL") == 0)
			{
				model_number = to!int(strip(line[5 .. $]));

				this.models[model_number] = new Model(model_number);

				model_new = true;
			}

			// Ignore non-atoms
			if (line.indexOf("ATOM") != 0)
				continue;

			// Check if a model is available. 
			if (!model_new)
			{
				this.models[model_number] = new Model(model_number);
				model_new = true;
			}

			// Filter for known residues
			auto residue_name = to!string(strip(line[17 .. 20]));
			string name = to!string(strip(line[12 .. 16])); // atom name
			if (!residue_codes.canFind(residue_name) || !atom_names.canFind(name))
				continue;

			// Check for new residue
			int curr_residue = to!int(strip(line[22 .. 26]));

			if (curr_residue != last_residue)
			{
				this.models[model_number].residues[curr_residue] = new Residue(curr_residue);
				this.models[model_number].residues[curr_residue].name = residue_name;

				last_residue = curr_residue;
			}

			// Process atom data
			this.models[model_number].residues[curr_residue].atoms[name] = new Atom(name);

			this.models[model_number].residues[curr_residue].atoms[name].serial = to!int(
				strip(line[6 .. 11]));

			this.models[model_number].residues[curr_residue].atoms[name].altLoc = line[16];
			this.models[model_number].residues[curr_residue].atoms[name].chainID = line[21];
			this.models[model_number].residues[curr_residue].atoms[name].iCode = line[26];

			this.models[model_number].residues[curr_residue].atoms[name].x = to!double(
				strip(line[30 .. 38]));
			this.models[model_number].residues[curr_residue].atoms[name].y = to!double(
				strip(line[38 .. 46]));
			this.models[model_number].residues[curr_residue].atoms[name].z = to!double(
				strip(line[46 .. 54]));

			this.models[model_number].residues[curr_residue].atoms[name].occupancy = to!double(
				strip(line[54 .. 60]));
			this.models[model_number].residues[curr_residue].atoms[name].tempFactor = to!double(
				strip(line[60 .. 66]));
		}
		bar.finish();

		double start = 1.0;
		double stop = to!double(this.models.length);

		int uCPUs = nCPUs;

		// For a single job for one model!!!
		if (abs(stop - start) == 0.0)
			uCPUs = 1;

		int[] a = new int[](uCPUs);
		int[] b = new int[](uCPUs);

		double step = (stop - start) / to!double(uCPUs);

		// Devide the jobs among the CPUs
		a[0] = to!int(start);
		b[0] = to!int(floor(start + step));
		foreach (i; 1 .. uCPUs)
		{
			a[i] = b[i - 1] + 1;
			b[i] = to!int(floor(start + step * to!double(i + 1)));
		}

		//writefln("PDB = %s", this.ID);
		auto jobs = new Range[](uCPUs);
		foreach (i; 0 .. uCPUs)
		{
			//writefln("Job No %2s calculating differential geometry for models from %4s to %4s.", i+1, a[i], b[i]);
			jobs[i] = Range(a[i], b[i]);
		}

		//stdout.flush();

		foreach (job; parallel(jobs))
		{
			job.xgeo(this.models);
		}
	}

	// Save a PDB file.
	void save_pdb(in string file_name)
	{
		File out_file = File(file_name, "w");
		auto model_keys = sort!("a < b")(this.models.keys);
		foreach (mk; model_keys)
		{
			out_file.writefln("MODEL %4d", mk);

			auto res_keys = sort!("a < b")(this.models[mk].residues.keys);
			foreach (rk; res_keys)
			{
				auto atom_keys = this.models[mk].residues[rk].atoms.keys;

				// Sort the atoms by their serial numbers
				string[int] atom_serial;
				foreach (ak; atom_keys)
				{
					atom_serial[this.models[mk].residues[rk].atoms[ak].serial] = ak;
				}
				auto atom_serial_keys = sort!("a < b")(atom_serial.keys);

				foreach (ask; atom_serial_keys)
				{
					auto ak = atom_serial[ask];

					out_file.writefln("ATOM  %5d %4s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f",
						this.models[mk].residues[rk].atoms[ak].serial,
						this.models[mk].residues[rk].atoms[ak].name,
						this.models[mk].residues[rk].name,
						this.models[mk].residues[rk].atoms[ak].chainID,
						this.models[mk].residues[rk].resSeq,
						this.models[mk].residues[rk].atoms[ak].x,
						this.models[mk].residues[rk].atoms[ak].y,
						this.models[mk].residues[rk].atoms[ak].z, 
						this.models[mk].residues[rk].atoms[ak].occupancy,
						this.models[mk].residues[rk].atoms[ak].tempFactor);
				}
			}
			out_file.writeln("TER");
		}
		out_file.writeln("END");
		out_file.close();
	}

	// Save a PDB file.
	void save_csv(in string file_name)
	{
		File out_file = File(file_name, "w");
		auto model_keys = sort!("a < b")(this.models.keys);
		foreach (mk; model_keys)
		{
			auto res_keys = sort!("a < b")(this.models[mk].residues.keys);
			foreach (rk; res_keys)
			{
				out_file.writefln("%03d,%3s,%1s,%04d,%7.4f,%7.4f",
					mk,
					this.models[mk].residues[rk].name,
					this.models[mk].residues[rk].atoms["CA"].chainID,
					this.models[mk].residues[rk].resSeq,
					this.models[mk].residues[rk].curvature,
					this.models[mk].residues[rk].torsion);
			}
		}
		out_file.close();
	}

}
