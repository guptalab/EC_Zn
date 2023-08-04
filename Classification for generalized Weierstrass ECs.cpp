/*
##########################################################################################
Title: Classification of generalized Weierstrass elliptic curve over Zn, n >= 2
Developers: Param Parekh, Paavan Parekh
Graduate Mentor: Sourav Deb
Mentor: Prof. Manish K Gupta
Version: 1.0
Website: https://www.guptalab.org/ecc
This is source code used for classification of generalized weierstrass elliptic curve over Zn.
It includes calculation of No. of non-singular curves, No. of classes of generalized 
weierstrass EC over Z_n, n >= 2.
########################################################################################## 
*/
#include <iostream>
#include <bits/stdc++.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
using namespace std;
#define ll long long int

ll gcd(ll a, ll b)
{
    if (b == 0)
        return a;
    return gcd(b, a % b);
}


string print_curve(int a1, int a2, int a3, int a4, int a6)
{
    string str = "";
    //str+= "y^2  + " + to_string(a1) + "xy + "+ to_string(a3) + "y = x^3 + " + to_string(a2) + "x^2 + " +to_string(a4) + "x + " + to_string(a6) + ",";
    string a1str="",a3str="",a2str="",a4str="",a6str="";
   
    if(a1) a1str="+ "+to_string(a1)+"xy ";
    if(a3) a3str="+ "+to_string(a3)+"y ";
    if(a2) a2str="+ "+to_string(a2)+"x<sup>2</sup> ";
    if(a4) a4str="+ "+to_string(a4)+"x ";
    if(a6) a6str="+ "+to_string(a6);
    str += "y<sup>2</sup>  " + a1str + a3str + "= x<sup>3</sup> "+ a2str + a4str + a6str;
   
    return str;
}

long long int discriminant(int a1, int a2, int a3, int a4, int a6, int p)
{
    long long int b2 = (((a1 % p) * a1 % p) % p + (4 * a2) % p) % p; // :)))) detected on 18/11/22
    long long int b4 = ((2 * a4) % p + ((a1 % p) * a3 % p) % p) % p;
    long long int b6 = ((a3 % p * a3 % p) % p + (4 * a6) % p) % p;
    long long int b8 = (((a1 % p) * a1 % p * a6) % p + (4 * (a2 % p) * (a6 % p)) % p + p - ((a1 % p) * (a3 % p) * a4) % p + ((a2 % p) * (a3 % p) * a3) % p + p - ((a4 % p) * a4) % p) % p;
    long long int delta = (3 * p - ((b2 % p) * (b2 % p) * (b8 % p)) % p - (8 * (b4 % p) * (b4 % p) * (b4 % p)) % p - (27 * (b6 % p) * b6 % p) % p + (9 * (b2 % p) * (b4 % p) * (b6 % p)) % p) % p;
    return delta;
}

map<pair<string, vector<int>>, int> all_curves_1(ll p)
{
    map<pair<string, vector<int>>, int> mp;
    int cnt = 0;
    for (int a1 = 0; a1 < p; a1++)
    {
        for (int a2 = 0; a2 < p; a2++)
        {
            for (int a3 = 0; a3 < p; a3++)
            {
                for (int a4 = 0; a4 < p; a4++)
                {
                    for (int a6 = 0; a6 < p; a6++)
                    {

                        if (gcd(discriminant(a1, a2, a3, a4, a6, p) % p, p) == 1)
                        {

                            cnt++;
                            vector<int> v;
                            v.push_back(a1);
                            v.push_back(a3);
                            v.push_back(a2);
                            v.push_back(a4);
                            v.push_back(a6);
                            mp[{print_curve(a1, a2, a3, a4, a6), v}] = 0;
                        }
                    }
                }
            }
        }
    }
    return mp;
}

ll modpow(ll x, ll n, ll p)
{
    if (n == 0)
        return 1 % p;
    ll y = 0;
    if (n % 2 == 0)
    {
        y = modpow(x, n / 2, p);
        y = (y * y) % p;
    }
    else
    {
        y = x % p;
        y = (y * modpow(x, n - 1, p)) % p;
    }
    return (ll)((y + p) % p);
}
ll x, y;
ll exgcd(ll a, ll b, ll &x, ll &y)
{
    if (b == 0)
    {
        x = 1;
        y = 0;
        return a;
    }
    ll x1, y1;
    ll d = exgcd(b, a % b, x1, y1);
    x = y1;
    y = x1 - y1 * (a / b);
    return d;
}
ll modinv(ll a, ll m)
{
    // return modpow(a,m-2,m);
    ll x, y;
    exgcd(a, m, x, y);
    return (x % m + m) % m;
}

pair<string, vector<int>> transformed_curve(int u, int r, int s, int t, int A, int B, int C, int D, int F, int p)
{
    ll u_ = modinv(u, p) % p;

    ll a1 = (((A % p + (2 * s) % p) % p) * u_ % p) % p;                                                                                                             // xy
    ll a3 = ((B % p + (r * A) % p + (2 * t) % p) % p * modpow(u_, 3, p)) % p;                                                                                       // y
    ll a2 = (((2 * p + C % p - (A * s) % p - (s * s) % p + (3 * r) % p) % p) % p * modpow(u_, 2, p) % p) % p;                                                       // x^2
    ll a4 = ((((2 * C * r) % p + D % p + 4 * p - (A * s * r) % p - (A * t) % p - (B * s) % p - (2 * s * t) % p + (3 * r * r) % p) % p) * modpow(u_, 4, p) % p) % p; // x
    ll a6 = (((F % p + (r * r * r) % p + (C * r * r) % p + 3 * p - (A * r * t) % p - (B * t) % p + (D * r) % p - (t * t) % p) % p) * modpow(u_, 6, p) % p) % p;     // constant
    vector<int> vec;
    vec.push_back(a1);
    vec.push_back(a3);
    vec.push_back(a2);
    vec.push_back(a4);
    vec.push_back(a6);

    if (gcd(discriminant(a1, a2, a3, a4, a6, p) % p, p) == 1)
    {
        string str = print_curve(a1,a2,a3,a4,a6);
        return {str, vec};
    }
    else
        return {"singular curve", {}};
}
ll phin(ll n)
{
    // calculate phi(n) function

    ll result = 0;
    for (int i = 1; i < n; i++)
    {
        if (gcd(i, n) == 1)
        {
            result++;
        }
    }
    return result;
}


pair<ll,vector<pair<ll,ll>>> numberofpoints(vector<int> v, ll p)
{
   ll a1 = v[0];
    ll a3 = v[1];
    ll a2 = v[2];
    ll a4 = v[3];
    ll a6 = v[4];
   

    ll count = 0;
    vector<pair<ll,ll>> points;
    for (ll x = 0; x < p; x++)
    {
        for(ll y = 0;y < p;y++)
        {
            ll left = ((y%p * y%p) % p + (a1*x*y%p + a3*y%p)%p)%p;
            ll right = (modpow(x,3,p) + (a2%p * modpow(x,2,p))%p + a4 * x + a6) % p;
            if(left == right)
            {count++;
            points.push_back({x,y});
            }
        }
    }
    return {count,points};
}

int main()
{
    ofstream file1;
   { 
    file1.open("ViewClassificationOfCurves.html");
    file1 << "<!DOCTYPE html>\n";
    file1 << "<html>\n";
    file1 << "<head>\n";
    file1 << "<title>Classification Of Elliptic Curves Over Zn</title>\n";
    file1 << "</head>\n";
    file1 << "<body>\n";
    file1 << "</body>\n";
    file1 << "<center>\n";
    file1 << "<table border=\"1\">\n";
    file1 << "<tr>\n";
    file1 << "<th>Zn</th>\n";
    file1 << "<th> #Non-Singular Curves </th>\n";
    file1 << "<th>#Classes</th>\n";
    file1 << "</tr>\n";
    file1 << "</center>\n";
    }

   

    for (ll p = 1; p <= 50; p++)
    {
        cout << "executing for "<< p << endl;
        file1 << "<tr>\n";

        if (mkdir("Zn_classification") == -1)
        cerr << "Error :  " << strerror(errno) << endl;

         ofstream file2;
        string str_file = "Z_"+to_string(p)+".html";
         file2.open("Zn_classification/"+str_file);
        { /////////
                file2 << "<!DOCTYPE html>\n";
                file2 << "<html>\n";
                file2 << "<head>\n";
                file2 << "<title>Classification Of Elliptic Curves Over Z"+to_string(p)+"</title>\n";
                file2 << "</head>\n";
                file2 << "<body>\n";
                file2 << "</body>\n";
                file2 << "<center>\n";
                file2 << "<table border=\"1\">\n";
                file2 << "<tr>\n";
                file2 << "<th>class no. </th>\n";
                file2 << "<th> # of curves </th>\n";
                file2 << "<th>class leader</th>\n";
                 file2 << "<th># Aut(E)</th>\n";
             file2 << "<th># points</th>\n";
                file2 << "</tr>\n";
                file2 << "</center>\n";

         ////////////
         }
        
        map<pair<string, vector<int>>, int> mp = all_curves_1(p);
        int class_no = 0;
        map<int,pair<int,string>> class_detail;
          map<int, vector<pair<string, vector<int>>>> class_curves;
        map<string, int> cleader_points;
        map<string, ll> cleader_aut;

        for (auto curve = mp.begin(); curve != mp.end(); curve++)
        {

            string ccurve = curve->first.first;
            vector<int> coeff = curve->first.second;
            int color = curve->second;

            if (color == 0)
            {
                class_no++;
                class_detail[class_no].second = ccurve;
                 cleader_points[ccurve] = numberofpoints(coeff,p).first;
                ll autcnt = 0;
                for (ll u = 0; u < p; u++)
                {
                    if (gcd(u, p) != 1)
                        continue;

                    for (ll s = 0; s < p; s++)
                    {
                        for (ll r = 0; r < p; r++)
                        {
                            for (ll t = 0; t < p; t++)
                            {

                                pair<string, vector<int>> tcurve = transformed_curve(u, r, s, t, coeff[0], coeff[1], coeff[2], coeff[3], coeff[4], p);

   if(tcurve.second == coeff)autcnt++;
                                if (tcurve.first == "singular curve")
                                    continue;
                                else
                                {
                                    if (mp[tcurve] == 0)
                                    {
                                        ll a1 = tcurve.second[0];
                                        ll a3 = tcurve.second[1];
                                        ll a2 = tcurve.second[2];
                                        ll a4 = tcurve.second[3];
                                        ll a6 = tcurve.second[4];

                                        mp[tcurve] = class_no;
                                        class_curves[class_no].push_back({tcurve.first, tcurve.second});
                                        class_detail[class_no].first++;
                                    }
                                }
                            }
                        }
                    }
                }
                 cleader_aut[ccurve] = autcnt;
            }
        }

     
        file1 << "<td style=\"text-align:center\"><a href=\"../classification/Zn_classification/"<<str_file<<"\">"<<p<<"</a></td>\n";
        file1 << "<td style=\"text-align:center\">" << mp.size() << "</td>\n";
        file1 << "<td style=\"text-align:center\">" << class_no << "</td>\n";
       
         string fldname = "classes_Z"+to_string(p);
                 const char* fld = fldname.c_str();
                 if (mkdir(fld) == -1)
                   cerr << "Error :  " << strerror(errno) << endl;

        
        for (auto it = class_detail.begin(); it != class_detail.end(); it++)
           {
            string class_file = "Z_"+to_string(p)+"_class_"+to_string(it->first)+".html";
             file2 << "</tr>\n"; 
             file2 << "<td style=\"text-align:center\"><a href=\"../"+fldname+"/"<<class_file<<"\">"<< it->first << "</a></td>\n";
             file2 << "<td style=\"text-align:center\">" << it->second.first << "</td>\n";
             file2 << "<td style=\"text-align:center\">" << it->second.second << "</td>\n";
                file2 << "<td style=\"text-align:center\">" << cleader_aut[it->second.second] << "</td>\n";
            file2 << "<td style=\"text-align:center\">" << cleader_points[it->second.second] << "</td>\n";
             file2 << "</tr>\n";
           }
           file2 << "</table>\n";
           file2 << "</html>\n";
           file2.close();

       
       
         for (int i = 1; i <= class_no; i++)
        {
               

                ofstream file3;
                string cls_file = "Z_"+to_string(p)+"_class_"+to_string(i)+".html";
                file3.open(fldname+"/"+cls_file);
                 { /////////
                file3 << "<!DOCTYPE html>\n";
                file3 << "<html>\n";
                file3 << "<head>\n";
                file3 << "<title>Elliptic Curves Over Z"+to_string(p)+" : class +"+to_string(i)+" </title>\n";
                file3 << "</head>\n";
                file3 << "<body>\n";
                file3 << "class leader: "<< class_detail[i].second <<"\n";
                file3 << "</body>\n";
                file3 << "<center>\n";
                file3 << "<table border=\"1\">\n";
                file3 << "<tr>\n";
                file3 << "<th>Sr no.</th>\n";
                file3 << "<th>class members</th>\n";
                   file3 << "<th># points</th>\n";
                file3 << "<th>points</th>\n";
                file3 << "</tr>\n";
                file3 << "</center>\n";

               ////////////
                  }
              int pp=1;
              for (auto ttcurve : class_curves[i])
                      {
                             file3 << "</tr>\n"; 
                             file3 <<  "<td style=\"text-align:center\">"<< pp << "</td>\n";
                              file3 << "<td style=\"text-align:center\">" << ttcurve.first << "</td>\n";
                pair<ll,vector<pair<ll,ll>>> ans = numberofpoints(ttcurve.second, p);
                file3 << "<td style=\"text-align:center\">" << ans.first << "</td>\n";
                 file3 << "<td style=\"text-align:center\">";
                for(auto x:ans.second)
                {
                    file3 << "(" << x.first << "," << x.second << ") ";
                }
                file3 << "</td>\n";
                file3 << "</tr>\n";
                             pp++;
                      }
              file3 << "</table>\n";
             file3 << "</html>\n";
              file3.close();
        }
        

        file1 << "</tr>\n";

    }
    file1 << "</table>\n";
    file1 << "</html>\n";
     file1.close();
    return 0;
}