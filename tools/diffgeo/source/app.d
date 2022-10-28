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
module main;

import pdb;

// Standard libs
import std.conv;
import std.stdio;
import std.getopt;
import std.parallelism;

void main(string[] args)
{
	string id;

	int cCPUS = 1;

	auto helpInformation = getopt(args,
		"id", "file id: id.pdb -> id.csv", &id,
		"threads", "number of threads", &cCPUS);

	if (helpInformation.helpWanted)
	{
		defaultGetoptPrinter("This program compute curvature and torsion for PDB files.",
			helpInformation.options);
	}

	int nCPUs = totalCPUs;
	if (cCPUS <= totalCPUs)
		nCPUs = cCPUS;

	writeln();
	writeln("DiffGeo Neural Network 0.1");
	writeln("Copyright (c) 2022 Rinaldo Wander Montalvão");
	writeln();
	writefln("Threads : %d", nCPUs);
	writefln("PDB id  : %s\n", id);

	try
	{
		auto pdb_file = new PDB(id, id ~ ".pdb", nCPUs);

		pdb_file.save_csv(id ~ ".csv");
	}
	catch (Exception e)
	{
		writeln();
		if (id.length == 0)
			writeln("Error: empty file id.");
		else
			writeln(e.msg);
	}
}
