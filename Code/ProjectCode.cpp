#include "consoledisplays.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <stack>
#include <string.h>

using namespace std;


template <typename dataType>
struct GraphData { // Directed Graph

	dataType source;
	dataType destination;
	int edgeWeight;

	void set(dataType A, dataType B, int d) {
		this->source = A;
		this->destination = B;
		this->edgeWeight = d;
	}

	GraphData<dataType>* resize(GraphData<dataType> *obj, int _size) {

		GraphData<dataType> *_temp = new GraphData<dataType>[_size + 1];
		for (int i = 0; i < _size; i++)
			_temp[i] = obj[i];

		delete obj;
		return _temp;
	}


};


template<typename dataT>
int sumEdgeWeights(GraphData<dataT> *arr, int _size) {

	int total_sum = 0;
	std::cout << "MST Edge Weight= ";
	for (int i = 0; i < _size; i++) {

		total_sum += arr[i].edgeWeight;
		if (i + 1 != _size) {
			std::cout << arr[i].edgeWeight << " + ";
		}
		else {
			std::cout << arr[i].edgeWeight;
		}

	}
	std::cout << " = " << total_sum << std::endl;
	return total_sum;
}



template<typename dataType>
void addVertex(std::vector<dataType> &Vertices, dataType newElement) {

	bool alreadyExists = false;

	for (int i = 0; i < Vertices.size(); i++) {
		if (Vertices[i] == newElement) {
			alreadyExists = true;
			break;
		}
	}
	if (alreadyExists == false)
		Vertices.push_back(newElement);

}

template<typename datatype>
void DisplayEdgesList(GraphData<datatype> *arr, int arr_s) {

	for (int i = 0; i < arr_s; i++) {
		std::cout << arr[i].source << " --- ";
		std::cout << arr[i].edgeWeight << " --- ";
		std::cout << arr[i].destination << std::endl;
	}

}

template<typename dataType>
void swapEdges(GraphData<dataType> &a, GraphData<dataType> &b) {
	GraphData<dataType> temp = a;
	a = b;
	b = temp;
	return;
}


template< typename dataType >
void sortEdgeList(GraphData<dataType> *arr, int arr_s) {

	for (int i = 0; i < arr_s - 1; i++) {

		for (int j = 0; j < arr_s - i - 1; j++) {

			if (arr[j].edgeWeight > arr[j + 1].edgeWeight)
				swapEdges(arr[j], arr[j + 1]);

		}

	}

}

template <typename dataType>
bool CheckDisjoint(GraphData<dataType> *Paired, int count, GraphData<dataType> obj) {

	if (Paired == NULL)
		return true;

	int index_source_occ = -1; // index of source occurence
	int index_source_x = -1;
	int index_source_y = -1;

	int index_dest_occ = -1; // index of source occurence
	int index_dest_x = -1;
	int index_dest_y = -1;

	for (int i = 0; i < count; i++) {
		if (obj.source == Paired[i].source) {
			index_source_occ = i;
			index_source_x = i;
			break;
		}
		else if (obj.source == Paired[i].destination) {
			index_source_occ = i;
			index_source_y = i;
			break;
		}
	}

	for (int i = 0; i < count; i++) {
		if (obj.destination == Paired[i].source) {
			index_dest_occ = i;
			index_dest_x = i;
			break;
		}
		else if (obj.destination == Paired[i].destination) {
			index_dest_occ = i;
			index_dest_y = i;
			break;
		}
	}


	if (index_dest_occ == -1 || index_source_occ == -1) {
		return true; // if either one of them is new vertex 	
	}
	else { // polygonal cycle can occur so need to check

   // Triangular cycle check here

		dataType key1 = NULL, key2 = NULL;

		if (index_source_x != -1 && index_dest_x != -1) {
			key1 = Paired[index_source_y].source;
			key2 = Paired[index_dest_y].source;
			if (Paired[index_source_y].source == Paired[index_dest_y].source) {
				return false;
			}
		}
		else if (index_dest_y != -1 && index_source_y != -1) {
			key1 = Paired[index_source_x].destination;
			key2 = Paired[index_dest_x].destination;
			if (Paired[index_source_x].destination == Paired[index_dest_x].destination) {
				return false;
			}
		}
		else if (index_source_x != -1 && index_dest_y != -1) {
			key1 = Paired[index_source_x].destination;
			key2 = Paired[index_dest_y].source;
			if (Paired[index_source_x].destination == Paired[index_dest_y].source) {
				return false;
			}
		}
		else if (index_source_y != -1 && index_dest_x != -1) {
			key1 = Paired[index_source_y].source;
			key2 = Paired[index_dest_x].destination;
			if (Paired[index_source_y].source == Paired[index_dest_x].destination) {
				return false;
			}
		}

		//else { // triangular cycle not found
			  // Now we check for Polygon Shape

		bool _polygonkey1 = false;
		bool _polygonkey2 = false;

		for (int i = 0; i < count; i++) {

			if (i == index_dest_occ || i == index_source_occ)
				continue;

			if (Paired[i].source == key1 || Paired[i].destination == key1) {
				_polygonkey1 = true;
			}
			if (Paired[i].source == key2 || Paired[i].destination == key2) {
				_polygonkey2 = true;
			}

		}
		if (_polygonkey1 == false && _polygonkey2 == false) {
			return true;
		}
		else if (_polygonkey1 == true && _polygonkey2 == true) {
			return false;
		}
		else {
			return true;
		}

	}


}



template <typename dataType>
bool checkVisitedPathMST(std::vector<dataType> path , dataType source) {

	for (int i = 0; i < path.size(); i++) {
		if (path[i] == source) {
			return true;
		}
	}
	return false;
}






template <typename dataType>
void Display_MSTPath(GraphData<dataType> *edges, int count_edges) {

	sortEdgeList(edges, count_edges);// sort the edges

	std::vector<dataType > path;

	dataType indegreeZero_S = edges[0].source;
	bool indegreezs = false;
	dataType  indegreeZero_D = edges[0].destination;



		path.push_back(indegreeZero_S);

		path.push_back(indegreeZero_D);
	
	while (path.size() != count_edges+1 )
	{
		for (int i = 1; i < count_edges; i++) {

			bool checkS = checkVisitedPathMST(path, edges[i].source);
			bool checkD = checkVisitedPathMST(path, edges[i].destination);
			if (checkS == true && checkD != true) {
				path.push_back(edges[i].destination);
			}
			else if (checkD == true && checkS != true) {
				path.push_back(edges[i].source);
			}

		}
	}

	std::cout << "\n\n\n MST Path : ";

	for (int i = 0; i < path.size(); i++) {
		std::cout << path[i] << "  ";
		if (i + 1 != path.size())
			std::cout << " -> ";
	}
	

	
}


template <typename dataType>
GraphData<dataType>* KruskalMST(GraphData<dataType> *arr, int arr_s, std::vector<dataType> &vertices, int num_vert) {

	sortEdgeList(arr, arr_s); // sorts all the edges List

	//DisplayEdgesList(arr,arr_s);

	int count = 0; //

	GraphData<dataType> *PairedVertices = NULL; // included in kruskal

	GraphData<dataType> *temp = NULL;//  = new GraphData<dataType>[num_vert - 1];

	for (int i = 0; i < arr_s; i++) {

		bool _ifDisjoint = CheckDisjoint(PairedVertices, count, arr[i]);

		if (_ifDisjoint == true && count < num_vert - 1) {

			PairedVertices = PairedVertices->resize(PairedVertices, count);
			PairedVertices[count].source = arr[i].source;
			PairedVertices[count].destination = arr[i].destination;
			PairedVertices[count].edgeWeight = arr[i].edgeWeight;

			temp = temp->resize(temp, count);
			temp[count].source = arr[i].source;
			temp[count].destination = arr[i].destination;
			temp[count].edgeWeight = arr[i].edgeWeight;
			++count;
		}

	}

	return temp;
}


template <typename dataType >
GraphData<dataType>* readGraphData
(std::string filename, GraphData<dataType> *_temp, int &_sizetemp, std::vector<dataType> &Vertices) {

	std::fstream _fileobject;
	_fileobject.open(filename.c_str(), std::ios::in);

	//GraphData<> *_temp = NULL; // int _sizetemp = 0;

	int captureV1 = INT_MAX, captureV2 = INT_MAX, captureEdge = INT_MAX;


	if (_fileobject.is_open()) {

		while (!_fileobject.eof())
		{
			_fileobject >> captureV1;
			_fileobject >> captureV2;
			_fileobject >> captureEdge;

			_temp = _temp->resize(_temp, _sizetemp);
			_temp[_sizetemp].source = captureV1;
			_temp[_sizetemp].destination = captureV2;
			_temp[_sizetemp].edgeWeight = captureEdge;

			addVertex(Vertices, captureV1); // adds a vertex to vertices vector if already not added
			addVertex(Vertices, captureV2);

			++_sizetemp;


		}

	}
	else {
		std::cout << "Failed to Open File " << filename << std::endl;
	}
	return _temp;
}


template<typename dataType >
void dispayGraphVertices(std::vector<dataType> vertices) {

	if (vertices.size() != 0) {
		for (int i = 0; i < vertices.size(); i++)
			std::cout << "[" << vertices[i] << "]\t";
		std::cout << std::endl;
	}

}

void KruskalAlgo(std::string filename) {


	int _edgelistsize = 0;

	GraphData<int > *EdgeList = NULL;
	std::vector<int > Vertices;

	EdgeList = readGraphData(filename.c_str(), EdgeList, _edgelistsize, Vertices); // Edges are loaded from file here




	if (EdgeList != NULL) {

		std::cout << "=---------<< Graph Data from File >>---------=" << std::endl;

		DisplayEdgesList(EdgeList, _edgelistsize);

		std::cout << "=---------<< -------------------- >>---------=" << std::endl;

		std::cout << "=---------<< Graph Vertices >>---------=" << std::endl;

		dispayGraphVertices(Vertices);

		std::cout << "=---------<< -------------------- >>---------=" << std::endl;

		GraphData<int > *_edgesKruskal = NULL;

		_edgesKruskal = KruskalMST(EdgeList, _edgelistsize, Vertices, Vertices.size());

		std::cout << "=---------<<   Kruskal's Edges  >>---------=" << std::endl;

		DisplayEdgesList(_edgesKruskal, Vertices.size() - 1);

		std::cout << "=---------<< -------------------- >>---------=" << std::endl;

		int KMST_Weight = sumEdgeWeights(_edgesKruskal, Vertices.size() - 1);

		Display_MSTPath(_edgesKruskal, Vertices.size() - 1);


	}





}


//------------------------------------------------------------------------------------=

template<typename dataType>
bool Check_ifJoined(int index, GraphData<dataType> *arr, GraphData<dataType> *PairedVertices, int count) {

	bool joined = false;

	dataType k1 = arr[index].source;
	dataType k2 = arr[index].destination;


	if (count != 0) {

		for (int i = 0; i < count; i++) {

			if (PairedVertices[i].source == k1 || PairedVertices[i].destination == k1) {
				return true;
			}
			else if (PairedVertices[i].source == k2 || PairedVertices[i].destination == k2) {
				return true;
			}


		}
		return false;

	}
	else {
		return true;
	}
}

bool Check_ifExcludedIndex(int i, std::vector<int > excludedIndex) {

	for (int j = 0; j < excludedIndex.size(); j++)
		if (i == excludedIndex[j])
			return true;

	return false;
}

template <typename dataT>
int minimumEdge(GraphData<dataT> *arr, int arr_s,
	GraphData<dataT> *PairedVertices, int count, std::vector<int > excludedIndex) {

	// note excludedIndex consists of all the indexes that are already taken
	std::vector<int > newExcluded = excludedIndex;


	int mini = INT_MAX;
	bool if_excluded = false;
	bool if_joined = false;
	bool if_disjoint = false;
	int index = -1;

	for (int i = 0; i < arr_s; i++) {

		// 1st : Check if index is excluded means Edge list excluding the previous connected pairs
		if_excluded = Check_ifExcludedIndex(i, newExcluded);

		if (if_excluded == false) {

			// 2nd: Check if the new pair can join with the previous connected pairs
			if_joined = Check_ifJoined(i, arr, PairedVertices, count);

			if (if_joined == true) {

				// 3rd : Check if new pair does not make cycle
				if_disjoint = CheckDisjoint(PairedVertices, count, arr[i]);

				if (if_disjoint == true) {

					//4: Check if its minimum also
					if (mini > arr[i].edgeWeight) {
						mini = arr[i].edgeWeight;
						index = i;
					}

				}


			}



		}
	}


	if (index == -1)
		throw "exception throw at catching minimum!";
	else
		return index;

}


template<typename dataType>
GraphData<dataType>* PrimMST(GraphData<dataType> *arr, int arr_s, std::vector<dataType> vertices, int num_vert) {

	int count = 0; //

	sortEdgeList(arr, arr_s);

	GraphData<dataType> *PairedVertices = NULL; // included in prim's 

	std::vector<int > _excludedIndex;

	for (int i = 0; i < vertices.size() - 1; i++) {

		int index = minimumEdge(arr, arr_s, PairedVertices, count, _excludedIndex);

		_excludedIndex.push_back(index);

		PairedVertices = PairedVertices->resize(PairedVertices, count);
		PairedVertices[count] = arr[index];

		++count;
	}


	return PairedVertices;

}



void PrimAlgo(std::string filename) {


	int _edgelistsize = 0;

	GraphData<int > *EdgeList = NULL;
	std::vector<int > Vertices;

	EdgeList = readGraphData(filename.c_str(), EdgeList, _edgelistsize, Vertices); // Edges are loaded from file here

	if (EdgeList != NULL) {

		std::cout << "=---------<< Graph Data from File >>---------=" << std::endl;

		DisplayEdgesList(EdgeList, _edgelistsize);

		std::cout << "=---------<< -------------------- >>---------=" << std::endl;

		std::cout << "=---------<< Graph Vertices >>---------=" << std::endl;

		dispayGraphVertices(Vertices);

		std::cout << "=---------<< -------------------- >>---------=" << std::endl;

		GraphData<int > *_edgesPrim = NULL;

		_edgesPrim = PrimMST(EdgeList, _edgelistsize, Vertices, Vertices.size());

		std::cout << "=---------<<  Prim's Edges  >>---------=" << std::endl;

		DisplayEdgesList(_edgesPrim, Vertices.size() - 1);

		std::cout << "=---------<< -------------------- >>---------=" << std::endl;

		int PMST_Weight = sumEdgeWeights(_edgesPrim, Vertices.size() - 1);

		Display_MSTPath(_edgesPrim,Vertices.size()-1);


	}

}



//---------------------------Bipartite or complete bipartite /////////////////////////



// The concept is same as python dicts 
template <typename T, typename S>
struct key_valueMap {
	T key;
	S value;
};

template<typename T>
struct ColorMap {

	key_valueMap <T, std::string > *_colorArr;
	int _colorarrsize;
	int _insertionP; // insertion pointer


	ColorMap(int _size) {
		this->_colorarrsize = _size;
		this->_colorArr = new key_valueMap <T, std::string >[_size];

		for (int i = 0; i < this->_colorarrsize; i++) {
			this->_colorArr[i].key = NULL;
			this->_colorArr[i].value = "";
		}

		this->_insertionP = 0;
	}

	bool add_valuecolor(T _key, std::string &color) {

		int index = -1;
		bool valueAlreadyExists = this->check_ifKeyValueExists(_key, index);

		if (valueAlreadyExists == false && index == -1) {
			this->_colorArr[this->_insertionP].key = _key;
			this->_colorArr[this->_insertionP].value = color;
			++this->_insertionP;
			return true;
		}
		else if (valueAlreadyExists == true && index != -1) {
			color = this->_colorArr[index].value;
			return false;
		}


	}


	std::string returnColor(T _key) {

		bool val = false; int index = -1;
		val = this->check_ifKeyValueExists(_key, index);

		if (index != -1 && val == true) {
			return this->_colorArr[index].value;
		}
		else {
			return "NULL";
		}


	}

	bool check_ifKeyValueExists(T _key, int &index) {

		bool val = false;
		// here we check if the key/vertex alreadys exists
		for (int i = 0; i < this->_colorarrsize; i++) {
			if (_key == this->_colorArr[i].key) {
				val = true;
				index = i;
			}
		}

		if (val == false && index == -1)
			return false;
		else
			return true;

	}

};


template <typename dataType>
class Bipartite_Graph {


private:

	std::vector<dataType > *AdjacenyList;
	std::vector<dataType > Vertices;
	GraphData< dataType > *arr;
	int arr_size;

public:


	//void set_DirectedGraph(GraphData<dataType> arr[], int size_arr, char vertices[], int num_vertices)
	//{


	//	AdjacenyList = new std::vector<char>[num_vertices];
	//	EdgeWeight = new std::vector<int >[num_vertices];

	//	char Source = '-';

	//	for (int i = 0; i < num_vertices; i++) {
	//		Source = vertices[i];
	//		Vertices.push_back(Source);
	//		for (int j = 0; j < size_arr; j++) {

	//			if (arr[j].source == Source || arr[j].destination == Source) {

	//				if (arr[j].source == Source)
	//					AdjacenyList[i].push_back(arr[j].destination);
	//				else
	//					AdjacenyList[i].push_back(arr[j].source);

	//				EdgeWeight[i].push_back(arr[j].edgeWeight);
	//			}


	//		}

	//	}

	//}

	Bipartite_Graph(std::string filename) {
		this->AdjacenyList = NULL;
		this->arr = NULL;
		this->arr_size = 0; // number of edges in graph

		if (filename != "") {

			this->readGraphData(filename);
			if (this->arr == NULL) {
				std::cout << "Data could not be fetched!" << std::endl;
			}
			else {
				this->set_UnDirectedGraph();
				bool status = this->check_ifBipartiteGraph();
			}
		}
		else {
			std::cout << "Empty File like " << filename << " cant be loaded" << std::endl;
		}

	}

	void set_UnDirectedGraph() {


		int num_vertices = this->arr_size;

		this->AdjacenyList = new std::vector<dataType>[num_vertices];

		char Source = '-';


		for (int i = 0; i < this->Vertices.size(); i++) {

			Source = this->Vertices[i];
			std::cout << "[ " << Source << "] \t";
			for (int j = 0; j < this->arr_size - 1; j++) {

				if (this->arr[j].source == Source) {

					std::cout << arr[j].destination << "->";
					AdjacenyList[i].push_back(arr[j].destination);

				}

			}
			std::cout << std::endl;
		}

	}

	void display() const {


		std::cout << "=----------<< Adjaceny List of Graph >>-----------=" << std::endl;


		for (int i = 0; i < this->Vertices.size(); i++) {
			std::cout << "[" << Vertices[i] << "] ->";
			for (int j = 0; j < this->AdjacenyList[i].size(); j++) {
				std::cout << "(" << this->AdjacenyList[i][j] << ")" << " ->";
			}
			std::cout << "NULL" << std::endl;
		}



	}


	bool check_ifBipartiteGraph() {

		ColorMap<dataType> colormapping(this->Vertices.size());
		//std::vector<dataType> visitedVertices;

		std::string colorInput = "Red";
		std::string colorOutput = "Blue";

		bool status = false; // , if_exists = false;
		//int index = -1;

		bool BiPartiteGraph = true;

		status = colormapping.add_valuecolor(Vertices[0], colorInput);


		for (int i = 0; i < this->Vertices.size(); i++) {

			if (BiPartiteGraph == true) {

				//std::cout << "Vertex Selected = " << Vertices[i] << "\tIts Color = " << colormapping.returnColor(Vertices[i]) << std::endl;
				//visitedVertices.push_back(Vertices[i]);
				colorInput = colormapping.returnColor(Vertices[i]);

				for (int j = 0; j < this->AdjacenyList[i].size(); j++) {

					if (colorInput == "Red")
						colorOutput = "Blue";
					else
						colorOutput = "Red";

					status = colormapping.add_valuecolor(this->AdjacenyList[i][j], colorOutput);

					if (status == false && colorOutput != colorInput) {
						//std::cout << "okie it was already colored opposite ..." << std::endl;
						//std::cout << "[" << this->AdjacenyList[i][j] << "] [ " << colormapping.returnColor(this->AdjacenyList[i][j]) << " ] " << std::endl;
					}
					else if (status == false && colorOutput == colorInput) {
						BiPartiteGraph = false;
						std::cout << "Graph is not Bipartite" << std::endl;
						//std::cout << "[" << this->AdjacenyList[i][j] << "] [ " << colormapping.returnColor(this->AdjacenyList[i][j]) << " ] " << std::endl;
						break;
					}
					else {
						//std::cout << "[" << this->AdjacenyList[i][j] << "] [ " << colormapping.returnColor(this->AdjacenyList[i][j]) << " ] " << std::endl;
					}

				}

			}
			else {
				break;
			}
		}

		if (BiPartiteGraph == false) {

			std::cout << "=------------<< ----------------------- >>-----------=" << std::endl;
			std::cout << "=------------<< Graph is not Bi Partite >>-----------=" << std::endl;
			std::cout << "=------------<< ----------------------- >>-----------=" << std::endl;

			return false;
		}
		else {


			std::vector<dataType> setA, setB; // setA = red,setB= blue 

			for (int i = 0; i < this->Vertices.size(); i++) {

				if (colormapping._colorArr[i].value == "Red") {
					setA.push_back(colormapping._colorArr[i].key);
				}
				else if (colormapping._colorArr[i].value == "Blue") {
					setB.push_back(colormapping._colorArr[i].key);
				}


			}


			std::cout << "=------------<< Bi Partite Vertices Set >>-----------=" << std::endl;
			displaySet(setA, " Set A ");
			displaySet(setB, " Set B ");


			bool is_completeBipartite = this->check_ifCompleteBiPartite(setA, setB);
			if (is_completeBipartite == true) {
				std::cout << "Graph is Complete Bi Partite in Nature" << std::endl;
			}

			std::cout << "=------------<< ----------------------- >>-----------=" << std::endl;
			return true;
		}





	}



private:


	bool check_ifCompleteBiPartite(std::vector<dataType> setA, std::vector<dataType> setB) {

		int TotalDegree = this->Vertices.size() - 1;

		int _checkforsetB = setA.size();
		int _checkforsetA = setB.size();


		dataType Source = '-';

		// 1st: Check if Elements of Set A Degree =no. of elements in set B
		for (int i = 0; i < setA.size(); i++) {

			Source = setA[i];

			for (int find = 0; find < this->Vertices.size(); find++) {

				if (Source == Vertices[find]) {

					if (this->AdjacenyList[find].size() != _checkforsetA || this->AdjacenyList[find].size() == 0) {
						return false;
					}
					else {
						break;
					}
				}

			}

		}

		// 2nd: Check if Elements of Set B Degree =no. of elements in set A
		for (int i = 0; i < setB.size(); i++) {

			Source = setB[i];

			for (int find = 0; find < this->Vertices.size(); find++) {

				if (Source == Vertices[find]) {

					if (this->AdjacenyList[find].size() != _checkforsetB || this->AdjacenyList[find].size() == 0) {
						return false;
					}
					else {
						break;
					}
				}

			}

		}

		return true; // if all members of set A has degree = no. of elements in set B and vice versa
		// then return true its a complete bi partite graph
	}

	void readGraphData(std::string filename) {

		std::fstream _fileobject;
		_fileobject.open(filename.c_str(), std::ios::in);

		dataType captureV1 = NULL, captureV2 = NULL;


		if (_fileobject.is_open()) {

			while (!_fileobject.eof())
			{
				_fileobject >> captureV1;
				_fileobject >> captureV2;

				this->arr = this->arr->resize(this->arr, this->arr_size);
				this->arr[this->arr_size].source = captureV1;
				this->arr[this->arr_size].destination = captureV2;

				addVertex(this->Vertices, captureV1); // adds a vertex to vertices vector if already not added
				addVertex(this->Vertices, captureV2);

				++this->arr_size;

				if (_fileobject.eof())
					break;

			}

		}
		else {
			std::cout << "Failed to Open File " << filename << std::endl;
		}

		return;
	}

	void displaySet(std::vector<dataType> setAB, std::string name) const {

		std::cout << std::endl;
		std::cout << name;
		std::cout << " = { ";
		for (int i = 0; i < setAB.size(); i++) {
			std::cout << setAB[i];
			if (i + 1 != setAB.size())
				std::cout << " , ";
		}
		std::cout << " } " << std::endl;

	}



};


//-------------------------------------------------------------------=


void Bipartite_GraphAlgo(std::string filename) {
	Bipartite_Graph<char> object(filename);
}
/*---------------------------Euler path And Circuit--------------------------
-----------------------------------------------------------------------------
-----------------------------------------------------------------------------
-----------------------------------------------------------------------------
-----------------------------------------------------------------------------
-----------------------------------------------------------------------------*/
class Euler {
	int edge;
	vector<int> *graph;
	int *inDegree;
	int *outDegree;

public:
	Euler() {

	}
	Euler(int x) {
		edge = x;
		graph = new vector<int>[x];
		inDegree = new int[x];
		outDegree = new int[x];

		for (int i = 0; i < x; i++) {
			inDegree[i] = 0;
			outDegree[i] = 0;
		}

	}
	//Fucntion to Create a Adjacency List Using Vectors.
	//inDegree Array hold the InDegree of edges of a Vertex
	//OutDegree Array hold the OutDegree of edges of a Vertex
	void addEdge(int from, int to) {
		graph[from].push_back(to);
		inDegree[to]++;
		outDegree[from]++;
	}
	//
	int getEularianPath() {

		if (!graphHasEularianPath())
			return NULL;
		cout << "Path is: ";
		printCircuit();
	}
	//This function checks if the Eularian Path Exits
	//If at most OutDegree - InDegree = 1 then it's a Eularian Path
	//If Indegree = Outdegree then it's a Eularian Circuit
	bool graphHasEularianPath() {
		if (edge == 0) {
			return false;
		}
		int startNode = 0, endNode = 0;
		for (int i = 0; i < edge; i++) {
			if (outDegree[i] - inDegree[i] > 1 || inDegree[i] - outDegree[i] > 1) {
				cout << "It's not a Eularian graph!" << endl;
				return false;
			}
			else if (outDegree[i] - inDegree[i] == 1) {
				cout << "Starting Point of Euler Path is: " << i << endl;
				startNode++;
			}
			else if (inDegree[i] - outDegree[i] == 1) {
				cout << "Ending Point of Euler Path is: " << i << endl;
				endNode++;
			}

		}
		if (endNode == 0 && startNode == 0) {
			cout << "It's a Eularian Circuit" << endl;
		}
		else if (endNode > 0 || startNode > 0) {
			cout << "It's a Eularian Path" << endl;

		}
		return (endNode == 0 && startNode == 0) || (endNode == 1 && startNode == 1);
	}
	int findingStartNode() {
		int start = 0;
		for (int i = 0; i < edge; i++) {
			// Unique starting node.
			if (outDegree[i] - inDegree[i] == 1)
				return i;
			// Start at a node with an outgoing edge.
			if (outDegree[i] > 0) {
				start = i;
			}

		}
		return start;
	}
	void printCircuit()
	{
		// adj represents the adjacency list of
		// the directed graph
		// edge_count represents the number of edges
		// emerging from a vertex
		unordered_map<int, int> edge_count;

		for (int i = 0; i < edge; i++)
		{
			//find the count of edges to keep track
			//of unused edges
			edge_count[i] = graph[i].size();
		}

		// Maintain a stack to keep vertices
		stack<int> curr_path;

		// vector to store final circuit
		vector<int> circuit;

		// start from any vertex
		curr_path.push(findingStartNode());
		int curr_v = findingStartNode(); // Current vertex

		while (!curr_path.empty())
		{
			// If there's remaining edge
			if (edge_count[curr_v])
			{
				// Push the vertex
				curr_path.push(curr_v);
				// Find the next vertex using an edge
				int next_v = graph[curr_v].back();
				// and remove that edge
				edge_count[curr_v]--;
				graph[curr_v].pop_back();
				// Move to next vertex
				curr_v = next_v;
			}

			// back-track to find remaining circuit
			else
			{
				circuit.push_back(curr_v);

				// Back-tracking
				curr_v = curr_path.top();
				curr_path.pop();
			}
		}

		// we've got the circuit, now print it in reverse
		for (int i = circuit.size() - 1; i >= 0; i--)
		{
			cout << circuit[i];
			if (i)
				cout << " -> ";
		}
	}
	//Printing the Eulairan Graph's Adjacency List,Indegrees and outDegrees.
	void print() {
		vector<int>::iterator it;
		for (int i = 0; i < edge; i++) {
			cout << "AdjList of Matrix " << i << ": ";
			for (it = graph[i].begin(); it != graph[i].end(); it++) {
				cout << *it << "->";
			}
			cout << endl;
		}
		for (int i = 0; i < edge; i++) {
			cout << "In and OutDegrees of " << i << endl;
			cout << inDegree[i] << " " << outDegree[i] << endl;
		}
	}

};
void EulerGraph(string filename) {
	ifstream file;
	file.open(filename);
	int x, y;
	bool check = 0;
	file >> x;
	file >> y;
	cout << x << " " << y << endl;
	Euler obj(x);
	if (file.is_open()) {
		while (!file.eof()) {
			file >> x;
			file >> y;
			if (x == y) {
				break;
			}
			else
				obj.addEdge(x, y);
		}
	}
	obj.print();
	obj.getEularianPath();
}
/*-----------------------------------------------------------------------------
-----------------------------------------------------------------------------
-----------------------------------------------------------------------------
-------------------------Euler Ends here--------------------------*/


template <typename dataType>
class Cycle_Graph {


private:

	std::vector<dataType > *AdjacenyList;
	std::vector<dataType > Vertices;
	GraphData< dataType > *arr; // all edges
	int arr_size;
	GraphData <dataType > *PairedVertices;
	int total_cycles;

public:


	Cycle_Graph(std::string filename) {
		this->AdjacenyList = NULL;
		this->arr = NULL;
		this->arr_size = 0; // number of edges in graph

		if (filename != "") {

			this->readGraphData(filename);
			this->set_DirectedGraph();
			this->display();
			


		}
		else {
			std::cout << "Empty File like " << filename << " cant be loaded" << std::endl;
		}

	}

	void set_DirectedGraph() {


		int num_vertices = this->arr_size;

		this->AdjacenyList = new std::vector<dataType>[num_vertices];

		char Source = '-';

		for (int i = 0; i < this->Vertices.size(); i++) {

			Source = this->Vertices[i];

			for (int j = 0; j < this->arr_size; j++) {

				if (this->arr[j].source == Source) {

					AdjacenyList[i].push_back(arr[j].destination);

				}

			}

		}

	}

	void display() const {


		std::cout << "=----------<< Adjaceny List of Graph >>-----------=" << std::endl;


		for (int i = 0; i < this->Vertices.size(); i++) {
			std::cout << "[" << Vertices[i] << "] ->";
			for (int j = 0; j < this->AdjacenyList[i].size(); j++) {
				std::cout << "(" << this->AdjacenyList[i][j] << ")" << " ->";
			}
			std::cout << "NULL" << std::endl;
		}



	}

	void count_Cycles() {



		int count_pv = 0; // count of paired vertices
		this->PairedVertices = NULL;
		//int count_paired = 0;
		this->total_cycles = 0;

		for (int i = 0; i < this->arr_size; i++) {

			bool found_Cycle = detectCycle(count_pv, this->arr[i]);


		}

		std::cout << "\nTotal Cycles in Graph: " << total_cycles << std::endl;



	}

	bool checkVisited(dataType source, std::vector<int > visited) {

		for (int i = 0; i < visited.size(); i++) {
			if (visited[i] == source) {
				return true;
			}
		}
		return false;
	}

private:

	bool detectCycle(int &count_pv, GraphData<dataType> object) {

		// lets find if the object.source already exists in the destination of some
		// already paired vertices
		int index_source = -1;
		index_source = this->check_alreadyExistsin_Destination(count_pv, object);

		if (index_source == -1) {
			// if it does not exist just add it in paired vertices and retun back
			this->PairedVertices = this->PairedVertices->resize(this->PairedVertices, count_pv);
			this->PairedVertices[count_pv].source = object.source;
			this->PairedVertices[count_pv].destination = object.destination;
			++count_pv;
			//std::cout << "Being Added: " <<object.source << "---" << object.destination << std::endl;
			return false;
		}
		else { // if yes the source exists then we try finding 

			bool mutuallyLinked = checkMutualLinked(index_source, object.destination);
			if (mutuallyLinked == true) {
				std::vector<dataType> selfcycle;
				selfcycle.push_back(object.source);
				selfcycle.push_back(object.destination);
				this->displayCycle(selfcycle);
				++this->total_cycles;
				return true;
			}


			// now the mission is to keep on finding the object.dest 
			// in paired vertices where it is a source to some destination
			dataType find_ = object.source;
			bool _findstatus = false;
			//std::cout << "To be found = "  <<  find_ << std::endl;
			for (int i = 0; i < count_pv; i++) {

				if (this->PairedVertices[i].source == object.destination) {
					std::vector<dataType> cyclicComponents;
					addVertex(cyclicComponents, object.destination);
					std::vector<bool> visited(count_pv, false);
					_findstatus = backTrack_Search(count_pv, this->PairedVertices[i].destination, object, cyclicComponents, visited);
					if (_findstatus == true) {
						this->displayCycle(cyclicComponents);
						++this->total_cycles;
					}
				}


			}


			this->PairedVertices = this->PairedVertices->resize(this->PairedVertices, count_pv);
			this->PairedVertices[count_pv].source = object.source;
			this->PairedVertices[count_pv].destination = object.destination;

			++count_pv;

			return _findstatus;
		}



	}

	bool backTrack_Search(int &count_pv, dataType hook, GraphData<dataType> search, std::vector<dataType> &cyclicComponents, std::vector<bool > &visited) {

		for (int i = 0; i < count_pv; i++) {
			if (this->PairedVertices[i].source == hook && visited[i] != true) {
				visited[i] = true;
				addVertex(cyclicComponents, hook);
				//std::cout << "Source = " << this->PairedVertices[i].source << std::endl;
				if (this->PairedVertices[i].destination == search.source) {
					addVertex(cyclicComponents, search.source);
					//++this->total_cycles;
					//std::cout << "Hooked Search Found = " << this->PairedVertices[i].destination << std::endl;
					return true;
				}
				else {

					if (check_alreadyExistin_inPath(cyclicComponents, this->PairedVertices[i].destination) == false) {
						return backTrack_Search(count_pv, this->PairedVertices[i].destination, search, cyclicComponents, visited);
					}
					else {
						cyclicComponents.clear();
						addVertex(cyclicComponents, search.destination);
						addVertex(cyclicComponents, this->PairedVertices[i].destination);
						return backTrack_Search(count_pv, this->PairedVertices[i].destination, search, cyclicComponents, visited);
					}
				}
			}
		}

		return false;
	}

	bool check_alreadyExistin_inPath(std::vector<dataType> path, dataType element) {

		for (int i = 0; i < path.size(); i++) {
			if (path[i] == element)
				return true;
		}
		return false;

	}

	int check_alreadyExistsin_Destination(int count_pv, GraphData<dataType> object) {

		dataType var = NULL;
		for (int i = 0; i < count_pv; i++) {
			var = this->PairedVertices[i].destination;
			if (var == object.source) {
				return i;
			}
		}
		return -1;
	}

	bool checkMutualLinked(int index, dataType dest) {

		if (this->PairedVertices[index].source == dest) {
			return true;
		}
		else {
			return false;
		}


	}



	void readGraphData(std::string filename) {

		std::fstream _fileobject;
		_fileobject.open(filename.c_str(), std::ios::in);

		dataType captureV1 = NULL, captureV2 = NULL;


		if (_fileobject.is_open()) {

			while (!_fileobject.eof())
			{
				_fileobject >> captureV1;
				_fileobject >> captureV2;

				this->arr = this->arr->resize(this->arr, this->arr_size);
				this->arr[this->arr_size].source = captureV1;
				this->arr[this->arr_size].destination = captureV2;

				addVertex(this->Vertices, captureV1); // adds a vertex to vertices vector if already not added
				addVertex(this->Vertices, captureV2);

				++this->arr_size;

				if (_fileobject.eof())
					break;

			}

		}
		else {
			std::cout << "Failed to Open File " << filename << std::endl;
		}

		return;
	}


	void displayCycle(std::vector<dataType > cycle) {

		if (cycle.size() != 0)
			std::cout << "\nCycle : ";

		for (int i = 0; i < cycle.size(); i++)
			std::cout << cycle[i] << "\t";

	}

};

void GraphCyclesAlgo(std::string filename) {

	Cycle_Graph<int> object(filename);
	object.count_Cycles();

}



/*-----------------------------------------------------------------------------
-----------------------------------------------------------------------------
-----------------------------------------------------------------------------
-------------------------Relations--------------------------*/
class Relation {
	int **arr;
	int *dom;
	int *coDom;
	int row;
	int sizeA;
	int sizeB;
public:
	Relation(int setA[], int setB[], int sizeA, int sizeB) {
		row = sizeA * sizeB;
		this->sizeA = sizeA;
		this->sizeB = sizeB;
		dom = new int[sizeA];
		coDom = new int[sizeB];
		arr = new int*[row];
		for (int i = 0; i < row; i++) {
			arr[i] = new int[2];
		}

		//Applying Cartesian Product to make the Ordered Pairs.
		int domain, codomain;
		domain = codomain = 0;
		for (int x = 0; x < row; x++) {
			if (codomain >= sizeB) {
				codomain = 0;
				domain++;
			}
			for (int y = 0; y < 2; y++) {
				if (y == 0) {
					arr[x][y] = setA[domain];
				}
				if (codomain < sizeB && y != 0) {

					arr[x][y] = setB[codomain];
					codomain++;
				}
			}

		}
		//displaying 2dArray
		//for (int x = 0; x < row; x++) {
		//	for (int y = 0; y < 2; y++) {
		//		cout << arr[x][y] << " ";
		//	}
		//	cout << endl;
		//}
		for (int i = 0; i < sizeA; i++) {
			dom[i] = setA[i];

		}
		for (int i = 0; i < sizeB; i++) {
			coDom[i] = setB[i];
		}

	}

	int numofRealtions() {

	}
	void isReflexive() //This Function Checks the Number of Relfexive Relation in the Set and Displays them
	{
		int count = 0;
		int check = INT_MAX;
		cout << "\n\n=-------------------<< Showing Reflexive Relations >>---------------=" << endl;
		for (int x = 0; x < row; x++) {
			check = INT_MAX;
			for (int y = 0; y < 2; y++) {
				if (check == arr[x][y]) {
					cout << "(" << check << "," << arr[x][y] << ")" << " is Reflexive" << endl;
					count++;
				}
				else
					check = arr[x][y];
			}

		}
		cout << "Number of Reflexive Relations are: ";
		cout << count;
		cout << endl;
		//cout << "-*-*-*-*-*-*-**-*-*-*-*-*-*Relfexive check Ends Here!-*-*-*-*-*-*-*-*-**-*-*-*-*";
		cout << endl;
	}
	void isSymmetric() { //This Function Checks the Number of Symmetric Relation in the Set and Displays them
		int count = 0;
		int i = INT_MAX;
		int j = INT_MAX;
		i = arr[0][0];
		j = arr[0][1];
		cout << "\n\n=-------------------<< Showing Symmetric Relations >>---------------=" << endl;
		for (int x = 0; x < row; x++) {
			i = arr[x][0];
			j = arr[x][1];
			count = count + checkSymmetric(i, j);
		}
		cout << "Number of Symmetric Relations Are :  ";
		cout << count << endl;
		//cout << "-*-*-*-*-*-*-**-*-*-*-*-*-*Symmetric check Ends Here!-*-*-*-*-*-*-*-*-**-*-*-*-*";
		cout << endl;
	}
	int checkSymmetric(int i, int j) { //An extended fucntion to Check the Symmetric Relation
		for (int y = 0; y < row; y++) {
			for (int z = 0; z < 2; z++) {
				if (z == 0) {
					if (i == arr[y][z + 1] && j == arr[y][z]) {
						cout << "(" << i << "," << j << ") " << "and (" << arr[y][z] << "," << arr[y][z + 1] << ")" << endl;
						return 1;
					}
				}

			}
		}
		return 0;

	}
	void isAntiSymmetric() {
		int count = 0;
		for (int i = 0; i < sizeA; i++) {
			for (int j = 0; j < sizeB; j++) {
				if (dom[i] == coDom[j]) {
					count++;
					if (count >= 1) {
						break;
					}
				}
			}
		}
		cout << "------------------------Showing Anti-Symmetric------------" << endl;
		if (count == 1) {
			for (int x = 0; x < row; x++) {
				int i = arr[x][0];
				int j = arr[x][1];
				count = count + checkSymmetric(i, j);
			}
		}
		else
			cout << "Anti-Symmetric Doesn't Hold" << endl;
	}


};
void RelationsAlgo() {
	int *setA, *setB;
	int sizeA = 0, sizeB = 0;
	cout << "\n\n\tEnter the number of elements in Set A: ";
	cin >> sizeA;
	cout << "\n\n\tEnter the number of elements in set B: ";
	cin >> sizeB;
	if (sizeA < 3 || sizeB < 3) {
		cout << "Number of elements should be greater than or equal to 3!" << endl;
		cout << "Enter again!" << endl;
		cout << "Enter the number of elements in Set A: " << endl;
		cin >> sizeA;
		cout << "Enter the number of elements in set B: " << endl;
		cin >> sizeB;
	}
	if (sizeA >= 3 || sizeB >= 3) {
		setA = new int[sizeA];
		setB = new int[sizeB];
		cout << "\n\n\tEnter the elements of Set A: " << endl;
		for (int i = 0; i < sizeA; i++) {
			std::cout << "\tA[" << i << "] : ";
			cin >> setA[i];
		}
		cout << "\n\n\tEnter the elements of Set B: " << endl;
		for (int i = 0; i < sizeB; i++) {
			std::cout << "\tB["<< i <<"] : ";
			cin >> setB[i];
		}
		Relation obj(setA, setB, sizeA, sizeB);
		obj.isReflexive();
		obj.isSymmetric();
		obj.isAntiSymmetric();
	}
}



int graph_menu() {

	int input;
	std::cout << "\n\t=----------<< Graphs >>------------=" << std::endl;
	std::cout << "\tPress 1==> PRIM MST" << std::endl;
	std::cout << "\tPress 2==> Kruskal MST" << std::endl;
	std::cout << "\tPress 3==> Bi Partite Graph" << std::endl;
	std::cout << "\tPress 4==> Euler Circuit/Path" << std::endl;
	std::cout << "\tPress 5==> Total cycle in the Graph" << std::endl;
	std::cout << "\tPress 6==> Exit" << std::endl;
	std::cout << "\n\t=----------<< ------ >>------------=" << std::endl;
	std::cout << "\n\tEnter respective no.: ";
	cin >> input;
	return input;

}

std::string enterthefilename() {
	std::string filename;
	std::cin.ignore();
	std::cout << "\n\n\tEnter the [File name.txt] to Read Data from: ";
	std::getline(std::cin, filename);
	std::cout << "\n\n" << std::endl;
	return filename;
}

void Graphs() {

	bool stop_ctrl = false;
	std::string filename;

	while (stop_ctrl != true)
	{
		system("Color 07");
		system("cls");
		int main_menu_value = graph_menu();
		switch (main_menu_value)
		{
		case 1:
			console_PRIMMST();
			filename = enterthefilename();
			PrimAlgo(filename);
			break;
		case 2:
			console_KRUSKALMST();
			filename = enterthefilename();
			KruskalAlgo(filename);
			break;
		case 3:
			console_BPGRAPH();
			filename = enterthefilename();
			Bipartite_GraphAlgo(filename);
			break;
		case 4:
			console_Euler();
			filename = enterthefilename();
			EulerGraph(filename);
			break;
		case 5:
			console_CylesGraph();
			filename = enterthefilename();
			GraphCyclesAlgo(filename);

		default:
			stop_ctrl = true;
			break;
		}



		system("pause");
	}
	//system("pause");
}

int main_menu() {

	int input;
	std::cout << "\n\t=----------<< Discrete Project >>------------=" << std::endl;
	std::cout << "\tPress 1==> Relations" << std::endl;
	std::cout << "\tPress 2==> Graphs" << std::endl;
	std::cout << "\tPress 3==> Credits" << std::endl;
	std::cout << "\tPress 4==> Exit" << std::endl;
	std::cout << "\n\t=----------<< ------ >>------------=" << std::endl;
	std::cout << "\n\tEnter respective no.: ";
	cin >> input;
	return input;

}

void Credits() {

	system("Color 0C");
	std::cout << "\n\n\n\t=----------<< Credits >>------------=\n" << std::endl;
	std::cout << "\t19F-0385 : M.Sameed Zahoor " << std::endl;
	std::cout << "\t19-0280  : Muhammad Ahmad " << std::endl;
	std::cout << "\t19F-0311 : Zaim Rana" << std::endl;
	std::cout << "\n\n\t=----------<< ----------- >>------------=" << std::endl;
}

int main(void) {

	bool stop_ctrl = false;
	std::string filename;

	while (stop_ctrl != true)
	{
		system("Color 07");
		system("cls");
		int main_menu_value = main_menu();
		switch (main_menu_value)
		{
		case 1:
			console_Relations();
			RelationsAlgo();
			break;
		case 2:
			Graphs();
			break;
		case 3:
			Credits();
			break;
		case 4:
			stop_ctrl = true;
			break;
		default:
			stop_ctrl = true;
			break;
		}

		system("pause");
	}

	system("pause");
}