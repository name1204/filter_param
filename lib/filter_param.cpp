/*
 * filter_param.cpp
 *
 *  Created on: 2021/07/22
 *      Author: matsu
 */

#include "filter_param.hpp"



using namespace std;

FILE* fileopen(const string& filename, const char mode, const string& call_file, const int call_line)
{
	FILE* fp = fopen(filename.c_str(), &mode);
	if(fp == NULL)
	{
		fprintf(stderr, "Error: [%s l.%d]Can't open file.(file name : %s, mode : %c)\n",
				call_file.c_str(), call_line, filename.c_str(), mode);
		exit(EXIT_FAILURE);
	}
	return fp;
}

vector<string> split_char(string& str, char spliter)
{
    stringstream ss{str};
    vector<std::string> strs;
    string buf;
    while (std::getline(ss, buf, spliter))
    {
        strs.emplace_back(buf);
    }
    return strs;
}

string BandParam::sprint()
{
	string str;
	switch(band_type)
	{
		case BandType::Pass:
		{
			auto fmt = string("PassBand(%6.3f, %6.3f)");
			str = format(fmt, left_side, right_side);
			break;
		}
		case BandType::Stop:
		{
			str = format("StopBand(%6.3f, %6.3f)", left_side, right_side);
			break;
		}
		case BandType::Transition:
		{
			str = format("TransitionBand(%6.3f, %6.3f)", left_side, right_side);
			break;
		}
	}
	return str;
}

FilterParam::FilterParam
(unsigned int pole, unsigned int zero, vector<BandParam> input_bands,
		unsigned int input_nsplit_approx, unsigned int input_nsplit_transition,
		double gd)
:n(pole), m(zero), bands(input_bands),
 nsplit_approx(input_nsplit_approx), nsplit_transition(input_nsplit_transition),
 group_delay(gd)
{
	// 周波数帯域の整合性チェック
	double band_left = 0.0;
	for(auto bp :bands)
	{
		if(bp.left() != band_left)
		{
			fprintf(stderr, "Error: [%s l.%d]Adjacent band edge is necessary same.\n"
					"Or first band of left side must be 0.0.\n",
					__FILE__, __LINE__);
			for(auto bp :bands)
			{
				fprintf(stderr, "%s\n", bp.sprint().c_str());
			}
			exit(EXIT_FAILURE);
		}
		band_left = bp.right();
	}
	if(bands.back().right() != 0.5)
	{
		fprintf(stderr, "Error: [%s l.%d]Last band of right side must be 0.5.\n",
				__FILE__, __LINE__);
		for(auto bp :bands)
		{
			fprintf(stderr, "%s\n", bp.sprint().c_str());
		}
		exit(EXIT_FAILURE);
	}
}

vector<FilterParam> FilterParam::read_csv(string& filename)
{
	vector<FilterParam> filter_params;
	ifstream ifs(filename);
	if(!ifs)
	{
		fprintf(stderr, "Error: [%s l.%d]Can't open file.(file name : %s, mode : read)\n",
				__FILE__, __LINE__, filename.c_str());
		exit(EXIT_FAILURE);
	}

	string buf;
	getline(ifs, buf);		// ヘッダ読み飛ばし
	while(getline(ifs, buf))
	{
		auto vals = split_char(buf, ',');
		auto type = FilterParam::analyze_type(vals.at(3));
		vector<BandParam> bands;

		switch(type)
		{
			case FilterType::LPF:
			{
				auto edges = FilterParam::analyze_edges(vals.at(3));
				if(edges.size() != 2)
				{
					fprintf(stderr, "Error: [%s l.%d]Format of filter state is illegal.(input : \"%s\")\n"
							"If you assign filter type L.P.F. , length of edges is only 2.\n",
							__FILE__, __LINE__, vals.at(3).c_str());
					exit(EXIT_FAILURE);
				}
				bands = FilterParam::gen_bands(FilterType::LPF, edges.at(0), edges.at(1));
				break;
			}
			default:
			{
				fprintf(stderr, "Error: [%s l.%d]It has not been implement yet.(input : \"%s\")\n",
						__FILE__, __LINE__, vals.at(3).c_str());
				exit(EXIT_FAILURE);
			}
		}

		filter_params.emplace_back
		(
			FilterParam
			(
				atoi(vals.at(1).c_str()), atoi(vals.at(2).c_str()),
				bands, atoi(vals.at(5).c_str()), atoi(vals.at(6).c_str()),
				atof(vals.at(4).c_str())
			)
		);
	}


	return filter_params;
}

/* フィルタのタイプを簡易入力する
 *   LPF(0.2, 0.3) : L.P.F.で通過域端 0.2, 阻止域端 0.3のフィルタ
 *
 *   Other(path) : その他の特性のフィルタ。pathは周波数帯域の情報を記したファイルへの相対パス
 */
FilterType FilterParam::analyze_type(string& input)
{
	FilterType type;
	string str = input;
	str.erase(remove(str.begin(), str.end(), ' '), str.end());
	string csv_space = regex_replace(str, regex("[(),:]+"), " ");
	auto input_type = split_char(csv_space, ' ');

	if(input_type.size() <= 1)
	{
		fprintf(stderr, "Error: [%s l.%d]Format of filter type is illegal.(input : \"%s\")\n",
				__FILE__, __LINE__, input.c_str());
		exit(EXIT_FAILURE);
	}
	else if(input_type.at(0) == "LPF")
	{
		type = FilterType::LPF;
	}
	else if(true)
	{
		fprintf(stderr, "Error: [%s l.%d]It has not been implement yet.(input : \"%s\")\n",
				__FILE__, __LINE__, input.c_str());
		exit(EXIT_FAILURE);
	}
	else
	{
		fprintf(stderr, "Error: [%s l.%d]Filter Type is undefined.(input : \"%s\")\n",
				__FILE__, __LINE__, input.c_str());
		exit(EXIT_FAILURE);
	}

	return type;
}

vector<double> FilterParam::analyze_edges(string& input)
{
	vector<double> edge;
	string str = input;
	str.erase(remove(str.begin(), str.end(), ' '), str.end());
	string csv_space = regex_replace(str, regex("[(),:]+"), " ");
	auto input_type = split_char(csv_space, ' ');
	input_type.erase(input_type.begin());

	if(input_type.size() <= 0)
	{
		fprintf(stderr, "Error: [%s l.%d]Format of filter state is illegal.(input : \"%s\")\n",
				__FILE__, __LINE__, input.c_str());
		exit(EXIT_FAILURE);
	}
	for(auto v :input_type)
	{
		edge.emplace_back(atof(v.c_str()));
	}
	return edge;
}
