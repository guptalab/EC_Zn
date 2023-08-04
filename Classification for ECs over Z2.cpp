/*
##########################################################################################
Title: Classification of Weierstrass elliptic curve over Z_2
Developers: Param Parekh, Paavan Parekh
Graduate Mentor: Sourav Deb
Mentor: Prof. Manish K Gupta
Version: 1.0
Website: https://www.guptalab.org/ecc
This is source code used for classification of elliptic curve over Z_2.
It includes calculation of No. of non-singular curves, No. of classes of EC over Z_2.
########################################################################################## 
*/

#include <iostream>
#include <bits/stdc++.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
using namespace std;
#define ll long long int
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
string print_curve_j_not0(ll a2, ll a6)
{
    string str = "";
    string a2str = "", a6str = "";
    if (a2)
        a2str = "+ " + to_string(a2) + "x<sup>2</sup> ";
    if (a6)
        a6str = "+ " + to_string(a6);

    str += "y<sup>2</sup> + xy = x<sup>3</sup> " + a2str + a6str;
    return str;
}
string print_curve_j0(ll a3, ll a4, ll a6)
{
    string str = "";
    string a3str = "", a4str = "", a6str = "";
    if (a3)
        a3str = "+ " + to_string(a3) + "y ";
    if (a4)
        a4str = "+ " + to_string(a4) + "x ";
    if (a6)
        a6str = "+ " + to_string(a6);

    str += "y<sup>2</sup> " + a3str + "= x<sup>3</sup> " + a4str + a6str;
    return str;
}
ll gcd(ll a, ll b)
{
    if (b == 0)
        return a;
    return gcd(b, a % b);
}
long long int discriminant_j_not0(ll a6, ll p)
{
    ll delta;
    delta = a6 % p;
    return delta % p;
}
long long int discriminant_j0(ll a3, ll p)
{
    ll delta;
    delta = modpow(a3, 4, p);
    return delta % p;
}
map<pair<string, vector<int>>, int> all_curves_1(ll p) // j!=0 char 2
{
    map<pair<string, vector<int>>, int> mp;
    for (ll a2 = 0; a2 < p; a2++)
    {
        for (ll a6 = 0; a6 < p; a6++)
        {
            if (gcd(discriminant_j_not0(a6, p) % p, p) == 1)
            {
                vector<int> v;
                v.push_back(a2);
                v.push_back(a6);
                mp[{print_curve_j_not0(a2, a6), v}] = 0;
            }
        }
    }
    return mp;
}
map<pair<string, vector<int>>, int> all_curves_2(ll p) // j==0 char 2
{
    map<pair<string, vector<int>>, int> mp;
    for (ll a3 = 0; a3 < p; a3++)
    {
        for (ll a4 = 0; a4 < p; a4++)
        {
            for (ll a6 = 0; a6 < p; a6++)
            {
                if (gcd(discriminant_j0(a3, p) % p, p) == 1)
                {
                    vector<int> v;
                    v.push_back(a3);
                    v.push_back(a4);
                    v.push_back(a6);
                    mp[{print_curve_j0(a3, a4, a6), v}] = 0;
                }
            }
        }
    }
    return mp;
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

pair<string, vector<int>> transformed_curve_j_not0(int s, int a2, int a6, int p)
{

    ll a2_ = ( p*p + a2 + (s % p * s % p) % p + s) % p;
    ll a6_ = a6;

    vector<int> vec;
    vec.push_back(a2_);
    vec.push_back(a6_);

    if (gcd(discriminant_j_not0(a6_, p) % p, p) == 1)
    {
        string str = print_curve_j_not0(a2_, a6_);
        return {str, vec};
    }
    else
        return {"singular curve", {}};
}
pair<string, vector<int>> transformed_curve_j0(int u, int s, int t, int a3, int a4, int a6, int p)
{
    int u_ = modinv(u, p);
    ll a3_ = (a3 % p * (modpow(u_, 3, p)) % p) % p;
    ll a4_ = ((modpow(s, 4, p) % p + (a3 % p * s % p) % p + a4 % p) % p * (modpow(u_, 4, p)) % p) % p;
    ll a6_ = ((t * t) % p + (a3 * t) % p + modpow(s, 6, p) % p + (a4%p * modpow(s,2,p)) % p + a6 % p) % p * (modpow(u_, 6, p)) % p; // put more mods :))

    vector<int> vec;
    vec.push_back(a3_);
    vec.push_back(a4_);
    vec.push_back(a6_);

    if (gcd(discriminant_j0(a3_, p) % p, p) == 1)
    {
        string str = print_curve_j0(a3_, a4_, a6_);
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
    ll a = v[0];
    ll b = v[1];
    ll count = 0;
    vector<pair<ll,ll>> points;
    for (ll x = 0; x < p; x++)
    {
        for(ll y = 0;y < p;y++)
        {
            if((y%p*y%p)%p == (modpow(x,3,p) + a*x + b)%p)
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
        file1.open("ViewClassificationOfCurvesChar2.html");
        file1 << "<!DOCTYPE html>\n";
        file1 << "<html>\n";
        file1 << "<head>\n";
        file1 << "<title>Classification Of Elliptic Curves Over Z<sub>2</sub></title>\n";
        file1 << "</head>\n";
        file1 << "<body>\n";
        file1 << "</body>\n";
        file1 << "<center>\n";
        file1 << "<h1>Classification Of Elliptic Curves Over Z<sub>2</sub>, j(E) != 0</h1>\n";
        file1 << "<table border=\"1\">\n";
        file1 << "<tr>\n";
        file1 << "<th>Zn</th>\n";
        file1 << "<th> #Non-Singular Curves </th>\n";
        file1 << "<th>#Classes</th>\n";
        file1 << "</tr>\n";
        file1 << "</center>\n";
    }
    for (ll p = 2; p <= 2; p++)
    {
        file1 << "<tr>\n";

        if (mkdir("Zn_classification") == -1)
            cerr << "Error :  " << strerror(errno) << endl;

        ofstream file2;
        string str_file = "Z_" + to_string(p) + "j_not0.html";
        file2.open("Zn_classification/" + str_file);
        { /////////
            file2 << "<!DOCTYPE html>\n";
            file2 << "<html>\n";
            file2 << "<head>\n";
            file2 << "<title>Classification Of Elliptic Curves Over Z<sub>" + to_string(p) + "</sub></title>\n";
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
        map<int, pair<int, string>> class_detail;
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
                for (ll s = 0; s < p; s++)
                {
                    pair<string, vector<int>> tcurve = transformed_curve_j_not0(s, coeff[0], coeff[1], p);

  if(tcurve.second == coeff)autcnt++;
                    if (tcurve.first == "singular curve")
                        continue;
                    else
                    {
                        if (mp[tcurve] == 0)
                        {
                            ll a2 = tcurve.second[0];
                            ll a6 = tcurve.second[1];

                            mp[tcurve] = class_no;
                             class_curves[class_no].push_back({tcurve.first, tcurve.second});
                            class_detail[class_no].first++;
                        }
                    }
                }
                cleader_aut[ccurve] = autcnt;
            }
        }

        file1 << "<td style=\"text-align:center\"><a href=\"../classification/Zn_classification/" << str_file << "\">" << p << "</a></td>\n";
        file1 << "<td style=\"text-align:center\">" << mp.size() << "</td>\n";
        file1 << "<td style=\"text-align:center\">" << class_no << "</td>\n";

        string fldname = "classes_Z" + to_string(p);
        const char *fld = fldname.c_str();
        if (mkdir(fld) == -1)
            cerr << "Error :  " << strerror(errno) << endl;

        for (auto it = class_detail.begin(); it != class_detail.end(); it++)
        {
            string class_file = "Z_" + to_string(p) + "_class_" + to_string(it->first) + "j_not0.html";
            file2 << "</tr>\n";
            file2 << "<td style=\"text-align:center\"><a href=\"../" + fldname + "/" << class_file << "\">" << it->first << "</a></td>\n";
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
            string cls_file = "Z_" + to_string(p) + "_class_" + to_string(i) + "j_not0.html";
            file3.open(fldname + "/" + cls_file);
            { /////////
                file3 << "<!DOCTYPE html>\n";
                file3 << "<html>\n";
                file3 << "<head>\n";
                file3 << "<title>Elliptic Curves Over Z" + to_string(p) + " : class +" + to_string(i) + " </title>\n";
                file3 << "</head>\n";
                file3 << "<body>\n";
                file3 << "class leader: " << class_detail[i].second << "\n";
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
            int pp = 1;
            for (auto ttcurve : class_curves[i])
            {
                file3 << "</tr>\n";
                file3 << "<td style=\"text-align:center\">" << pp << "</td>\n";
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




     /////////////////////////////////////////////////////////////////
     /////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////
    file1 << "<center>\n";
    file1 << "<h1>Classification Of Elliptic Curves Over Z<sub>2</sub>, j(E) = 0</h1>\n";
    file1 << "<table border=\"1\">\n";
    file1 << "<tr>\n";
    file1 << "<th>Zn</th>\n";
    file1 << "<th> #Non-Singular Curves </th>\n";
    file1 << "<th>#Classes</th>\n";
    file1 << "</tr>\n";
    file1 << "</center>\n";
    ///////////////////////////////////////////////
    for (ll p = 2; p <= 2; p++)
    {
        file1 << "<tr>\n";

        if (mkdir("Zn_classification") == -1)
            cerr << "Error :  " << strerror(errno) << endl;

        ofstream file2;
        string str_file = "Z_" + to_string(p) + "j0.html";
        file2.open("Zn_classification/" + str_file);
        { /////////
            file2 << "<!DOCTYPE html>\n";
            file2 << "<html>\n";
            file2 << "<head>\n";
            file2 << "<title>Classification Of Elliptic Curves Over Z<sub>" + to_string(p) + "</sub></title>\n";
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

        map<pair<string, vector<int>>, int> mp = all_curves_2(p);
        int class_no = 0;
        map<int, pair<int, string>> class_detail;
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
                        for (ll t = 0; t < p; t++)
                        {

                            pair<string, vector<int>> tcurve = transformed_curve_j0(u, s, t, coeff[0], coeff[1], coeff[2], p);

 if(tcurve.second == coeff)autcnt++;
                            if (tcurve.first == "singular curve")
                                continue;
                            else
                            {
                                if (mp[tcurve] == 0)
                                {
                                    mp[tcurve] = class_no;
                                      class_curves[class_no].push_back({tcurve.first, tcurve.second});
                                    class_detail[class_no].first++;
                                }
                            }
                        }
                    }
                }
                cleader_aut[ccurve] = autcnt;
            }
        }

            file1 << "<td style=\"text-align:center\"><a href=\"../classification/Zn_classification/" << str_file << "\">" << p << "</a></td>\n";
            file1 << "<td style=\"text-align:center\">" << mp.size() << "</td>\n";
            file1 << "<td style=\"text-align:center\">" << class_no << "</td>\n";

            string fldname = "classes_Z" + to_string(p);
            const char *fld = fldname.c_str();
            if (mkdir(fld) == -1)
                cerr << "Error :  " << strerror(errno) << endl;

            for (auto it = class_detail.begin(); it != class_detail.end(); it++)
            {
                string class_file = "Z_" + to_string(p) + "_class_" + to_string(it->first) + "j0.html";
                file2 << "</tr>\n";
                file2 << "<td style=\"text-align:center\"><a href=\"../" + fldname + "/" << class_file << "\">" << it->first << "</a></td>\n";
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
                string cls_file = "Z_" + to_string(p) + "_class_" + to_string(i) + "j0.html";
                file3.open(fldname + "/" + cls_file);
                { /////////
                    file3 << "<!DOCTYPE html>\n";
                    file3 << "<html>\n";
                    file3 << "<head>\n";
                    file3 << "<title>Elliptic Curves Over Z" + to_string(p) + " : class +" + to_string(i) + " </title>\n";
                    file3 << "</head>\n";
                    file3 << "<body>\n";
                    file3 << "class leader: " << class_detail[i].second << "\n";
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
                int pp = 1;
                for (auto ttcurve : class_curves[i])
                {
                    file3 << "</tr>\n";
                    file3 << "<td style=\"text-align:center\">" << pp << "</td>\n";
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

        file1 << "</html>\n";
        file1.close();
        return 0;
    }