/*
###############################################################################################
Title: Data analysis for conjecture : "The number of non singular reduced Weierstrass 
       elliptic curves over Z_(p^m), gcd(p^m,6)=1 is Phi(p^(2m))"
Developers: Param Parekh, Paavan Parekh
Graduate Mentor: Sourav Deb
Mentor: Prof. Manish K Gupta
Version: 1.0
Website: https://www.guptalab.org/ecc
This is source code used for the analysis of conjecture : "The number of non singular 
reduced Weierstrass elliptic curves over Z_(p^m), gcd(p^m,6)=1 is Phi(p^(2m))". Reduced 
weierstrass Elliptic curve y^2 = x^3 + ax + b is non-singular iff discriminant (Delta) 
is unit over Zpm, gcd(p^m,6)=1. So to find number of non-singular elliptic curves over 
Z_(p^m), we need to find how many pairs of (a,b) gives discriminant as unit element
Another way to look at this is as follows and this approach is what we have implemented:
#ways discriminant is unit = (total curves over Z_(p^m)) - (#ways discriminant is non-unit) 
                           = p^(2m) - (#ways discriminant is non-unit) 
where discriminant = 16(4a^3+27b^2),  a,b in Z_(p^m).
############################################################################################### 
*/

#include <bits/stdc++.h>
#define ll long long int
using namespace std;
ll gcdExtended(ll a, ll b, ll *x, ll *y)
{

    // Base Case
    if (a == 0)
    {
        *x = 0, *y = 1;
        return b;
    }

    // To store results of recursive call
    ll x1, y1;
    int gcd = gcdExtended(b % a, a, &x1, &y1);

    // Update x and y using results of recursive
    // call
    *x = y1 - (b / a) * x1;
    *y = x1;

    return gcd;
}
ll modpow(ll x, ll y, ll p)
{
    ll res = 1;
    while (y > 0)
    {
        if (y % 2 == 1)
            res = (res * x);
        y = y >> 1;
        x = (x * x);
    }
    return res % p;
}
ll modInverse(ll A, ll M)
{
    ll x, y;
    ll g = gcdExtended(A, M, &x, &y);
    if (g != 1)
        return -1;
    else
    {
        ll res = (x % M + M) % M;
        return res;
    }
}
ll Phi(ll n)
{
    ll res = 0;
    for (ll i = 0; i < n; i++)
    {
        if (__gcd(i, n) == 1)
            res++;
    }
    return res;
}

int main()
{
    
    fstream file1;
    ll p;
    cout << "Enter the value of p: ";
    cin >> p;
    for (int m = 1; m <= 4; m++)
    {
        ll N = modpow(p, m, LONG_MAX);
        //make folder p_data using mkdir command
        system(("mkdir "+to_string(p)+"_data").c_str());

        //open file in write mode
        file1.open(to_string(p)+"_data\\delta_analysis_" + to_string(N) + ".html", ios::out);
        
        file1 << "<!DOCTYPE html>\n";
        file1 << "<html>\n";
        file1 << "<head>\n";
        file1 << "<title>Data Analysis for Number of Solution to delta = non-unit</title>\n";
        file1 << "</head>\n";
        file1 << "<body>\n";
        file1 << "</body>\n";
        file1 << "<center><b> N = " << to_string(N) << "</b></center>\n";
        file1 << "<center>\n";
        file1 << "<table border=\"1\">\n";
        file1 << "<tr>\n";
        file1 << "<th>[C1] non-units</th>\n";
        file1 << "<th>[C2] #non-units lead to same <br>#solution for <img src=\"https://latex.codecogs.com/svg.image?\\Delta&space;\" alt=\"discriminant\"></th>\n";
        file1 << "<th>[C3] #solution for <img src=\"https://latex.codecogs.com/svg.image?\\Delta&space;\" alt=\"discriminant\"><br> = non-units in C1</th>\n";
        
        file1 << "</tr>\n";
        file1 << "</center>\n";

        unordered_map<ll, ll> analyser;
        for (ll a = 0; a < N; a++)
        {
            for (ll b = 0; b < N; b++)
            {
                ll delta = ((N * N - 16) % N * ((4 % N * modpow(a, 3, N) % N) % N + (27 % N * modpow(b, 2, N) % N) % N) % N) % N;
                // cout << delta << a << b << endl;
                if (__gcd(delta, N) != 1)
                    analyser[delta]++;
            }
        }
        unordered_map<ll, pair<ll, vector<ll>>> analyser2;

        for (auto i : analyser)
        {
            analyser2[i.second].first++;
            analyser2[i.second].second.push_back(i.first);
        }
        ll sum = 0;
        for (auto it : analyser2)
            sum += (it.first) * (it.second.second.size());

        for (auto i : analyser2)
        {
            file1 << "<tr>\n";
            file1 << "<td style=\"text-align:center\">";
            // loop till second last element via iterator, not auto
            for (auto j = i.second.second.begin(); j != i.second.second.end() - 1; j++)
            {
                file1 << to_string(*j) << ", ";
            }
            // print last element
            if (i.second.second.size() != 0)
                file1 << to_string(*(i.second.second.end() - 1));

            file1 << "</td>" << endl;
            file1 << "<td style=\"text-align:center\">" << to_string(i.second.second.size()) << "</td>" << endl;
            file1 << "<td style=\"text-align:center\">" << to_string(i.first) << "</td>" << endl;
            file1 << "</tr>\n";
        }
        
        file1 << "<tr>\n";
        file1 << "<td style=\"text-align:center\"></td>" << endl;
        file1 << "<td style=\"text-align:center;background-color:#fad263\"> C<sub>2</sub><sup>T</sup>C<sub>3</sub> = #ways delta can be non-unit</td>" << endl;
        file1 << "<td style=\"text-align:center;background-color:#fad263\">" << to_string(sum) << "</td>" << endl;
        file1 << "</tr>\n";
        file1 << "</table>\n";
        file1 << "Total non-singular EC over Z<sub>" << to_string(N) << "</sub> = " << to_string(N)<<"<sup>2</sup> - #ways delta can be non-unit (C<sub>2</sub><sup>T</sup>C<sub>3</sub>) = "<< to_string(modpow(N,2,LONG_MAX) - sum) << endl;
        file1 << "</html>\n";
        file1.close();
    }
}
