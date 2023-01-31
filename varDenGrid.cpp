#include <bits/stdc++.h>

using namespace std;

vector<vector<float>> arr;
vector<vector<float>> dummy_arr;
vector<vector<float>> scaledArr;
vector<vector<float>> distances;
vector<int> cluster;
vector<int> final_cluster;
vector<string> grid;
set<string> hypercubes;
int cluster_number = 1;

// Function to calculate distance between 2 points
float distance(vector<float> arr1, vector<float> arr2)
{
    float dist = 0;
    for(int i=0; i<arr1.size(); i++)
        dist = dist + pow((arr1[i] - arr2[i]) , 2);
    dist = sqrt(dist);
    return dist;
}

float algorithm(int no_of_clusters, int grids, string s){
    for(int i=0; i<dummy_arr.size(); i++){
        for(int j=0; j<dummy_arr[0].size(); j++)
            arr[i][j] = dummy_arr[i][j];
    }
    for(int i=0; i<cluster.size(); i++)
    cluster[i]=-1;
    cluster_number=1;
    hypercubes.clear();
    string output_distance_filename = "distance_" + s;
    string output_hypercube_filenmae = "hypercube_" + s;
    ofstream outfile(output_distance_filename);
    ofstream gridfile(output_hypercube_filenmae);

    // Missing Value Evaluation and scaling
    for(int i=0; i<arr[0].size() ; i++)
    {
        float max=0 , min=100000;

        // Calculating max and min
        for(int k=0 ; k<arr.size() ; k++)
        {
            if(arr[k][i]>max)
            max = arr[k][i];

            if(arr[k][i]<min && arr[k][i]>=0)
            min = arr[k][i];
        }
        int count = 0;
        // Missing Value Evaluation
        for(int p=0 ; p<arr.size() ; p++)
        {
            if(arr[p][i] == -1)
            {
            arr[p][i] = max+1;
            count++;
            }
        }
        if(count > 0)
        {
        max = max + 1;
        }
        // Scaling
        for(int j=0 ; j<arr.size() ; j++)
        {
            scaledArr[j][i] = ( arr[j][i] - min ) / ( max - min) * grids;
        }
    }

    // Making the distance matrix
    for(int i=0; i<distances.size(); i++){
        for(int j=0; j<distances.size(); j++){
            distances[i][j] = distance(scaledArr[i], scaledArr[j]);
        }
    }

    for(int lp=0; lp<no_of_clusters; lp++){

    // Assigning hypercubes
    for(int i=0; i<scaledArr.size(); i++){
        string s = "";
        if(cluster[i] == -1){
        for(int j=0; j<scaledArr[0].size(); j++){
            int temp = floor(scaledArr[i][j]);
            if(scaledArr[i][j]==(double)temp && scaledArr[i][j] != 0)
            temp = temp - 1;

            s = s+ to_string(temp);
        }
        grid[i] = s;
        hypercubes.insert(s);
        }
    }

    int max = 0;
    string max_cube;


    // Writing hypercubes in file and output
    for(auto itr : hypercubes){
        gridfile << "Hypercube: " << itr << "\n";
        // cout << "Hypercube: " << itr << "\n";
        int c = 0;
        for(int i=0; i<scaledArr.size(); i++){
        if(cluster[i] == -1){
            if(itr == grid[i]){
                c++;
                for(int j=0; j<scaledArr[0].size(); j++){
                    // cout << scaledArr[i][j] << " ";
                    gridfile << scaledArr[i][j] << " ";
                }
                // cout << "\n";
                gridfile << " Hypercube: " << grid[i] << "\n";
            }
        }
        }
        if(c>max) // Max hypercubes
        {
        max = c;
        max_cube = itr;
        }
        // cout << "\n";
        gridfile << "\n";
    }


    vector<int> hypercube_points;

    // Saving which points are in hypercube
    int x = 0;
    for(int i=0; i<scaledArr.size(); i++){
        if(cluster[i] == -1 && max_cube == grid[i]){
            hypercube_points.push_back(i);
            cluster[hypercube_points[x]] = cluster_number;
            x++;
        }
    }

    // Calculating density of max hypercube
    float e = 0;
    for(auto itr : hypercubes){
        if(itr == max_cube){
            for(int i=0; i<max; i++){
                float min = 1e10;
                for(int j=0; j<max; j++){
                    int x = hypercube_points[i];
                    int y = hypercube_points[j];
                    if(distances[x][y]<min && i!=j)
                    {
                    min = distances[x][y];
                    }
                // cout << "Distances: " << distances[x][y] << " ";
                }
                // cout << "Min: " << min << "\n";
                if(min > e)
                {
                e = min;
                // cout << "E: " << e << "\n";
                }
            }
            
            if(e == 1e10)
            e = 0;
            
            break;
        }
    }

    // Copying points in a vector for comparison
    vector<int> new_hypercube_points;

    for(int i=0; i<hypercube_points.size(); i++){
        new_hypercube_points.push_back(hypercube_points[i]);
    }

    // Forming a cluster
    vector<int> new_new_hypercube_points;
    int len;
    do{
        try{
        for(int i=0; i<scaledArr.size(); i++){
            if(cluster[i] == -1){
                double dist = 10000000;
                for(int j=0; j<new_hypercube_points.size(); j++){
                    if(distances[new_hypercube_points[j]][i] <= dist)
                    dist = distances[new_hypercube_points[j]][i];
                }
                if(dist <= e){
                new_new_hypercube_points.push_back(i);
                }
            }
        }
        }
        catch(bad_alloc & ba){
            break;
        }
        new_hypercube_points.clear();
        for(int i=0; i<new_new_hypercube_points.size(); i++){
            hypercube_points.push_back(new_new_hypercube_points[i]);
            new_hypercube_points.push_back(new_new_hypercube_points[i]);
            cluster[new_new_hypercube_points[i]] = cluster_number;
        }
        len = new_new_hypercube_points.size();
        new_new_hypercube_points.clear();
    } while(len > 0);

    cout << "Total non-empty hypercubes: " << hypercubes.size() << "\n";
    gridfile << "Total non-empty hypercubes: " << hypercubes.size() << "\n";
    cout << "Grid with max hypercube : " << max_cube << " with " << max << " points\n";
    gridfile << "Grid with max hypercube : " << max_cube << " with " << max << " points\n";
    cout << "Density of hypercube is : " << e << "\n";
    gridfile << "Density of hypercube is : " << e << "\n";

    // for(int i=0; i<hypercube_points.size(); i++){
    //     for(int j=0; j<scaledArr[0].size(); j++){
    //     cout << scaledArr[hypercube_points[i]][j] << " ";
    //     }
    //     cout << "\n";
    // }

    cout << "Number of points in cluster " << cluster_number << ": " << hypercube_points.size() << "\n";
    gridfile << "Number of points in cluster " << cluster_number << ": " << hypercube_points.size() << "\n";
    cluster_number ++;
    grid.clear();
    hypercubes.clear();
    hypercube_points.clear();
    new_hypercube_points.clear();
    new_new_hypercube_points.clear();

    }

    for(int i=0; i<scaledArr.size(); i++){
        if(cluster[i] == -1){
            double minn = 10000000;
            int cn=-1;
            double dist = 100000000;
            for(int j=0; j<scaledArr.size(); j++){

                if(cluster[j] != -1){
                dist = distances[i][j];
                if(dist < minn){
                minn = dist;
                cn = cluster[j];
                }
                }
            }
            cluster[i] = cn;
        }
    }
    outfile.close();
    gridfile.close();

    float intra_dist=0;
    float inter_dist=1e7;

    for(int i=0; i<scaledArr.size(); i++){
        for(int j=0; j<scaledArr.size(); j++){
            if(cluster[i]!=cluster[j]){
                if(distances[i][j]<inter_dist)
                inter_dist=distances[i][j];
            }
        }
    }
    for(int i=0; i<scaledArr.size(); i++){
        for(int j=0; j<scaledArr.size(); j++){
            if(cluster[i]==cluster[j] && i!=j){
                if(distances[i][j]>intra_dist)
                intra_dist=distances[i][j];
            }
        }
    }
    cout << "Inter: " << inter_dist << " " << "Intra: " << intra_dist << "\n";
    float dunn_index=inter_dist/intra_dist;

    return dunn_index;

    // for(int c=1; c<=no_of_clusters; c++){
    //     vector<float> centroid(scaledArr[0].size());
    // for(int i=0; i<scaledArr.size(); i++){
    //     if(cluster[i]==c){
            
    //     }
    // }
    // }
}

int main(int argc, char **argv){
    string s;
    s = argv[1];
    s.append(".txt");

    string mytext,temp = "";

    ifstream file(s);

    // Input the dataset
    while(getline(file, mytext)) {
        mytext = mytext + ' ';
        vector<float> row;
        for(int i=0; i<mytext.length(); i++) {
            if(mytext[i]!=' ') {
                temp = temp + mytext[i];
            }
            else
            {
                if(temp == "?")
                {
                temp = "-1";
                }
                row.push_back(stof(temp));
                temp="";
            }
        }
        dummy_arr.push_back(row);
    }
    // int o=0,zz=0;

    // for(int i=0; i<dummy_arr.size(); i++){
    //     if(dummy_arr[i][dummy_arr[0].size()-1] == 1)
    //     o++;
    //     else
    //     zz++;
    // }

    // cout << "One: " << o << "\n";
    // cout << "Two: " << zz << "\n";

    // Remove last column
    if(s=="wisconsin.txt" || s=="hepatitis.txt"){
        for(int i=0; i<dummy_arr.size(); i++){
            dummy_arr[i].pop_back();
        }
    }

    int class_num = 8;


    // Allocating space for the scaled vector
    scaledArr.resize(dummy_arr.size());
    for(int i=0; i<dummy_arr.size(); i++)
    {
    scaledArr[i].resize(dummy_arr[0].size());
    }

    // Allocating space for the dummy vector
    arr.resize(dummy_arr.size());
    for(int i=0; i<dummy_arr.size(); i++)
    {
    arr[i].resize(dummy_arr[0].size());
    }

    // Allocating space for the distances matrix
    distances.resize(dummy_arr.size());
    for(int i=0; i<dummy_arr.size(); i++)
    {
    distances[i].resize(dummy_arr.size());
    }

    // Allocating space to grid
    grid.resize(dummy_arr.size());

    // Allocating space to cluster
    cluster.resize(dummy_arr.size(), -1);

    // Allocating space to final cluster
    final_cluster.resize(dummy_arr.size(), -1);


    // Copying values in dummy arr
    for(int i=0; i<dummy_arr.size(); i++){
        for(int j=0; j<dummy_arr[0].size(); j++)
            arr[i][j] = dummy_arr[i][j];
    }

    int z=3;

    float dunn_index1 = algorithm(class_num,2,s);

    for(int i=0; i<cluster.size(); i++)
    final_cluster[i] = cluster[i];

    float dunn_index2 = algorithm(class_num,3,s);

    while(dunn_index2 >= dunn_index1 && z<=10){
    dunn_index1 = dunn_index2;

    for(int i=0; i<cluster.size(); i++)
    final_cluster[i] = cluster[i];

    z++;

    dunn_index2 = algorithm(class_num,z,s);
    }

    cout << "Value of z is: " << z-1 << " and Dunn index: " << dunn_index1 << "\n";

    int one=0,two=0,three=0,four=0,five=0,six=0,seven=0,eight=0;

    for(int i=0; i<cluster.size(); i++)
    {
        if(cluster[i]==1)
        one++;
        else if(cluster[i]==2)
        two++;
        else if(cluster[i]==3)
        three++;
        else if(cluster[i]==4)
        four++;
        else if(cluster[i]==5)
        five++;
        else if(cluster[i]==6)
        six++;
        else if(cluster[i]==7)
        seven++;
        else
        eight++;
    }

    cout << "One: " << one << "\n";
    cout << "Two: " << two << "\n";
    cout << "Three: " << three << "\n";
    cout << "four: " << four << "\n";
    cout << "five: " << five << "\n";
    cout << "six: " << six << "\n";
    cout << "seven: " << seven << "\n";
    cout << "eight: " << eight << "\n";
    
    file.close();

    return 0;
}