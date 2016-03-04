#include "smile/smile.h"
#include <cassert>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <limits>
#include <map>
#include <vector>
#include <iomanip>

using namespace std;

template < typename T > std::string to_string( const T& n ) {
  std::ostringstream stm ;
  stm << n ;
  return stm.str() ;
}

void printCPT(DSL_node *node)
{
  DSL_network* net = node->Network(); // node network                                                                   
  int handle = node->Handle();
  DSL_nodeDefinition *def = node->Definition();
  const DSL_Dmatrix &cpt = *def->GetMatrix();
  const DSL_idArray &outcomes = *def->GetOutcomesNames();
  const DSL_intArray &parents = net->GetParents(handle);
  int parentCount = parents.NumItems();
  
  DSL_intArray coords;
  for (int elemIdx = 0; elemIdx < cpt.GetSize(); elemIdx ++)
  {
    cpt.IndexToCoordinates(elemIdx, coords);
    cout << "P(" << node->GetId() << " = " << outcomes[coords[parentCount]] << " | ";
    for (int parentIdx = 0; parentIdx < parentCount; parentIdx ++)
    {
      if (parentIdx > 0) cout << ", ";
      DSL_node *parentNode = net->GetNode(parents[parentIdx]);
      const DSL_idArray &parentStates = *parentNode->Definition()->GetOutcomesNames();
      cout << parentNode->GetId() << " = " << parentStates[coords[parentIdx]];             
    }
    cout << ") = " << cpt[elemIdx] << endl;
  }
}

//CreateNetwork(children, bins, pos_probs, neg_probs, argv[3]);
void CreateNetwork( map< string, vector<string> > children, 
    map< string, vector<string> > bins,
    map< string, vector<string> > pos_probs,
    map< string, vector<string> > neg_probs,
    const char* xdsl_filename ) {

    int ret;
    DSL_network theNet;
    DSL_idArray outcomes;
//    DSL_stringArray outcomes;
    int uhandle, ohandle;
    // define node for each uid
    for( map< string, vector<string> >::iterator it = bins.begin(); it != bins.end(); ++it){
      // unobserved nodes
      // create node
      uhandle = theNet.AddNode(DSL_CPT,((it->first)).c_str());
      // set number of outcomes
      outcomes.Add("yes");
      ret = outcomes.Add("no");
      theNet.GetNode(uhandle)->Definition()->SetNumberOfOutcomes(outcomes);
      outcomes.Flush();
//      outcomes.Delete(1);//hacky
//      outcomes.Delete(2);//hacky
      // fill probability for leaf nodes (nodes without children)
      if( children[(it->first)].empty() ) {
        DSL_doubleArray leafProbs;
        leafProbs.SetSize( 2 );
        leafProbs[0] = 0.25;
        leafProbs[1] = 0.75;
        //leafProbs[0] = 0.0001;
        //leafProbs[1] = 0.9999;
        //leafProbs[0] = 0.5;
        //leafProbs[1] = 0.5;
        theNet.GetNode(uhandle)->Definition()->SetDefinition(leafProbs);
      }
      // observed nodes
      
      // create node
      ohandle = theNet.AddNode(DSL_CPT,((it->first)+"O").c_str());

      // set number of outcomes

      /*
      for( vector<string>::iterator it2 = (it->second).begin(); it2 != (it->second).end(); ++it2){
        //string name = "p_"+to_string(cnt);
        //string name = "p_"+to_string( round(atof((*it2).c_str())*10000)); // hacky..
        string name = (*it2);
        cout << name << endl;
        ret = outcomes.Add(name.c_str());
        cout << ret << endl;
      }*/
//      cout << (it->second).size() << endl;
      for(unsigned int index = 0 ; index < (it->second).size() ; ++index){
        string name = "bin" + to_string(index);
        ret = outcomes.Add(name.c_str());
        //cout << ret << endl;
        assert(ret == 0);
      }

//      printCPT(theNet.GetNode(ohandle));
      //cout << outcomes[1] << endl;
      theNet.GetNode(ohandle)->Definition()->SetNumberOfOutcomes(outcomes);

//      printCPT(theNet.GetNode(ohandle));

      outcomes.Flush();

      //create arc from unobserved node to observed node
      theNet.AddArc(uhandle, ohandle);

      // now fill in the conditional distribution for observed nodes
      // probability values from pos_probs and neg_probs
      DSL_doubleArray theProbs;
      theProbs.SetSize( (it->second).size() * 2 );
      for( unsigned int index = 0; index < (it->second).size() ; ++index ){
        if ( pos_probs.find((it->first)) == pos_probs.end() ) {
            // not found
            throw 20;
            exit(1);
        }

        try{
        theProbs[index] = atof((pos_probs[ (it->first) ].at(index)).c_str());
        theProbs[index + (it->second).size()] = atof((neg_probs[ (it->first) ].at(index)).c_str());
        }catch(int e){
          cerr << "dag and evidence doesn't match: " << it->first << endl;
          exit(1);
        }
      }
      theNet.GetNode(ohandle)->Definition()->SetDefinition(theProbs);
    }
cout << "blah" << endl;
    // create arcs (or edges) network (or ontology) structure
    int from_handle, to_handle;
    for( map< string, vector<string> >::iterator it = children.begin(); it != children.end(); ++it){  
      to_handle = theNet.FindNode((it->first).c_str());
      //TODO: check if to_handle is error code instead of actual handle!!! -2 if not found;;; these checkes should be done early on .. 
      if(to_handle == -2){

        // check if brenda head node 'whole body'
        if( (it->first).compare("BTO0001489") == 0 ){
          to_handle = theNet.AddNode(DSL_CPT,((it->first)).c_str());
          // set number of outcomes
          outcomes.Add("yes");
          ret = outcomes.Add("no");
          theNet.GetNode(to_handle)->Definition()->SetNumberOfOutcomes(outcomes);
          outcomes.Flush();
        }else{
	
          //cerr << (it->first) << " not found in evidence file" << endl;
          //exit(1);
          continue;
        }
      }

      for( vector<string>::iterator it2 = (it->second).begin(); it2 != (it->second).end(); ++it2){
        from_handle = theNet.FindNode((*it2).c_str());
        theNet.AddArc(from_handle, to_handle);
      }

      // now file in the conditional distribution for unobserved nodes
      // probability values from network structure
      if(! (it->second).empty() ) {
      
 //       printCPT(theNet.GetNode(to_handle));
        DSL_doubleArray theProbs;
        
        //cout << (it->second).size() << endl;
        DSL_Dmatrix *theCpt;
        theNet.GetNode(to_handle)->Definition()->GetDefinition(&theCpt);
        //theProbs.SetSize( ((it->second).size() + 1) * 2 );
        theProbs.SetSize( theCpt->GetSize() );
        for( unsigned int index = 0; index < theCpt->GetSize() - 2; index += 2){
          theProbs[index] = 1;
          theProbs[index+1] = 0;
        }
        theProbs[theCpt->GetSize() - 2] = 0.25; //
        theProbs[theCpt->GetSize() - 1] = 0.75; //

        //theProbs[theCpt->GetSize() - 2] = 0.0001; //
        //theProbs[theCpt->GetSize() - 1] = 0.9999; //
        //theProbs[theCpt->GetSize() - 2] = 0.0; //
        //theProbs[theCpt->GetSize() - 1] = 1.0; //


        theNet.GetNode(to_handle)->Definition()->SetDefinition(theProbs);
//        printCPT(theNet.GetNode(to_handle));
      }else{
//        printCPT(theNet.GetNode(to_handle));
      }
    }

    cout << "print network" << endl;

    theNet.WriteFile(xdsl_filename);
  };

void load_DAG(const char* DAG_filename, map<string, vector<string> >& children) {
   ifstream ifile(DAG_filename);

   string n1, n2;
   while(ifile >> n1 >> n2) {
       children[n1].push_back(n2);
   }
   ifile.close();
}

std::vector<string> &mysplit(const std::string &s, char delim, std::vector<string> &elems) {
    std::stringstream ss(s);
    string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}
std::vector<string> mysplit(const std::string &s, char delim) {
    std::vector<string> elems;
    mysplit(s, delim, elems);
    return elems;
}

void load_evidence(const char* evidence_filename, map<string, vector<string> >& bins, map<string, vector<string> >& pos_probs, map<string, vector<string> >& neg_probs) {
   ifstream ifile(evidence_filename);

   string uid;
   string b, pos_str, neg_str;
   while(ifile >> uid >> b >> pos_str >> neg_str) {
       bins[uid] = mysplit(b, '|');
       pos_probs[uid] = mysplit(pos_str, '|');
       neg_probs[uid] = mysplit(neg_str, '|');
   }
   ifile.close();
}

void myprint(map<string, vector<string> >& mymap){
    for( map<string, vector<string> >::iterator it = mymap.begin(); it != mymap.end(); ++it){
      cout << it->first << " ";
      for( vector<string>::iterator it2 = (it->second).begin(); it2 != (it->second).end(); ++it2){
        cout << *it2 << "|";
      }
      cout << endl;
    }
}

int main(int argc, char* argv[]) {

     if (argc < 4) {
       cout << "usage: train <DAG_filename> <evidence_filename> <xdsl_output_filename>" << endl;
       return -1;
     }

     map<string, vector<string> > children;
     map<string, vector<string> > bins, pos_probs, neg_probs;
cout << "load input" << endl;
     load_DAG(argv[1], children);
     load_evidence(argv[2], bins, pos_probs, neg_probs);
//myprint(pos_probs);
cout << "construct network" << endl;
     CreateNetwork(children, bins, pos_probs, neg_probs, argv[3]);

/*     
     map<int, int> go2node;
     DSL_network net;
     create_network(net, go2node, children, prob, prob_obs, argv[4]);
     
     map<int, double> result = prob;
     for(map<int, double>::iterator it = result.begin(); it != result.end(); ++it) {
       it->second = -1;
     }
     inference_with_network(net, go2node, result);
     
     ofstream ofile(argv[3]);
     for(map<int, double>::iterator it = result.begin(); it != result.end(); ++it) {
       ofile << setfill('0') << setw(7) << it->first << "\t" << it->second << endl;
     }
     ofile.close();
*/
}

