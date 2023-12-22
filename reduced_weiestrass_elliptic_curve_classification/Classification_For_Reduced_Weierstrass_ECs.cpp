/*
##########################################################################################
Title: Classification of reduced Weierstrass elliptic curve over Zn, gcd(n,6)=1, n >= 5
Developers: Param Parekh, Paavan Parekh
Graduate Mentor: Sourav Deb
Mentor: Prof. Manish K Gupta
Version: 1.0
visit https://www.guptalab.org/ecc for more information on results and research work
This is source code used for classification of reduced weierstrass elliptic curve over Zn.
It includes calculation of No. of non-singular curves, No. of classes of reduced 
weierstrass EC over Z_n, gcd(n,6)=1, n >= 5.
########################################################################################## 
*/


#include <iostream>
#include <bits/stdc++.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
using namespace std;
#define ll long long int

bool isprime(ll n)
{
    if (n == 1)
        return false;
    if (n == 2)
        return true;
    if (n % 2 == 0)
        return false;
    for (ll i = 3; i * i <= n; i += 2)
    {
        if (n % i == 0)
            return false;
    }
    return true;
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
string print_curve(ll a, ll b)
{
    string str = "";
    string astr = "", bstr = "";

    if (a)
        astr = "+ " + to_string(a) + "x ";
    if (b)
        bstr = "+ " + to_string(b);

    str += "y<sup>2</sup> = x<sup>3</sup> " + astr + bstr;

    return str;
}
ll gcd(ll a, ll b)
{
    if (b == 0)
        return a;
    return gcd(b, a % b);
}
long long int discriminant(ll a, ll b, ll p)
{
    ll delta;
    delta = ((p * p - 16) * ( (4%p * modpow(a,3,p)%p)%p + (27 * (b % p) * (b % p))%p ) % p) % p;
    return delta%p;
}
map<pair<string, vector<int>>, int> all_curves_1(ll p)
{
    map<pair<string, vector<int>>, int> mp;
    
    for (ll a = 0; a < p; a++)
    {
        for (ll b = 0; b < p; b++)
        {
            if (gcd(discriminant(a, b, p) % p, p) == 1)
            {
                vector<int> v;
                v.push_back(a);
                v.push_back(b);
                mp[{print_curve(a, b), v}] = 0;
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
pair<string, vector<int>> transformed_curve(int u, int A, int B, int p)
{
    ll u_ = modinv(u, p) % p;

    ll a = (A % p * modpow(u_, 4, p)) % p;
    ll b = (B % p * modpow(u_, 6, p)) % p;

    vector<int> vec;
    vec.push_back(a);
    vec.push_back(b);

    if (gcd(discriminant(a, b, p) % p, p) == 1)
    {
        string str = print_curve(a, b);
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
string get_current_dir_name()
{
    char buff[FILENAME_MAX];
    getcwd(buff, FILENAME_MAX);
    string current_working_dir(buff);
    string::size_type pos = current_working_dir.find_last_of("\\/");
    return current_working_dir.substr(pos + 1);
}
int main()
{
    ofstream file1;
    {
        file1.open("ViewClassificationOfReducedEllipticCurves.html");
        file1 << "<!DOCTYPE html>\n";
        file1 << "<html>\n";
        file1 << "<head>\n";
        file1 << "<title>Classification Of Elliptic Curves Over Z<sub>n</sub></title>\n";
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
    for (ll p = 4; p <= 30; p++)
    {
        if(p%2 == 0 || p%3 == 0 )continue;

        file1 << "<tr>\n";

        if (mkdir("Zn_classification") == -1)
        {
             if(errno != EEXIST)
            {
                cerr << "Error :  " << strerror(errno) << endl;
                return 0;
            }
        }

        ofstream file2;
        string str_file = "Z_" + to_string(p) + ".html";
        file2.open("Zn_classification/" + str_file);
        { /////////
            file2 << "<!DOCTYPE html>\n";
            file2 << "<html>\n";
            file2 << "<head>\n";
            file2 << "<title>Classification Of Elliptic Curves Over Z" + to_string(p) + "</title>\n";
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
                for (ll u = 0; u < p; u++)
                {
                    if (gcd(u, p) != 1)
                        continue;

                    pair<string, vector<int>> tcurve = transformed_curve(u, coeff[0], coeff[1], p);
                    
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
                cleader_aut[ccurve] = autcnt;
            }
        }
        string currdir = get_current_dir_name();
        file1 << "<td style=\"text-align:center\"><a href=\"../"<<currdir<<"/Zn_classification/" << str_file << "\">" << p << "</a></td>\n";
        file1 << "<td style=\"text-align:center\">" << mp.size() << "</td>\n";
        file1 << "<td style=\"text-align:center\">" << class_no << "</td>\n";

        string fldname = "classes_Z" + to_string(p);
        const char *fld = fldname.c_str();
        if (mkdir(fld) == -1)
        {
             if(errno != EEXIST)
            {
                cerr << "Error :  " << strerror(errno) << endl;
                return 0;
            }
        }

        for (auto it = class_detail.begin(); it != class_detail.end(); it++)
        {
            string class_file = "Z_" + to_string(p) + "_class_" + to_string(it->first) + ".html";
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
            string cls_file = "Z_" + to_string(p) + "_class_" + to_string(i) + ".html";
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
    file1 << "</html>\n";
    file1.close();
    return 0;
}