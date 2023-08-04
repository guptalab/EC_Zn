/*
#############################################################################################
Title: Data Analysis for conjectures "The number of quadratic and cubic residues along with 
       their occurrences over Z_(p^m)"
Developers: Param Parekh, Paavan Parekh
Graduate Mentor: Sourav Deb
Mentor: Prof. Manish K Gupta
Version: 1.0
Website: https://www.guptalab.org/ecc
This is source code used for Data Analysis for conjectures "The number of quadratic and cubic 
        residues along with their occurrences over Z_(p^m)". Simple bruteforce approach is 
        used to analyse the nature of both quadratic and cubic residues. Resultant data can
        be used to answer following questions:

        > How many quadratic residues are there over Z_(p^m) ?
        > How many cubic residues are there over Z_(p^m) ?
        > What can be said on the occurrence of quadratic residues ? 
        > What can be said on the occurrence of cubic residues ? 
        > How many quadratic residues are also cubic residues or How many six-residues are 
          there ? 
        > What can be said on the occurrence of six-residues (6Rs) ?   
############################################################################################# 
*/


#include <bits/stdc++.h>
using namespace std;
#define ll long long int

ll mod_pow(ll a, ll b, ll m)
{
    ll res = 1;
    while (b > 0)
    {
        if (b & 1)
            res = (res * a) % m;
        a = (a * a) % m;
        b >>= 1;
    }
    return res;
}

ll mod_inverse(ll a, ll m)
{
    ll m0 = m;
    ll y = 0, x = 1;

    if (m == 1)
        return 0;

    while (a > 1)
    {
        ll q = a / m;
        ll t = m;
        m = a % m, a = t;
        t = y;
        y = x - q * y;
        x = t;
    }
    if (x < 0)
        x += m0;
    return x;
}

ll power(ll a, ll b)
{
    ll res = 1;
    while (b > 0)
    {
        if (b & 1)
            res = (res * a);
        a = (a * a);
        b >>= 1;
    }
    return res;
}

int main()
{
    ll p, mp;
    cout << "Enter p: ";
    cin >> p;

    string dir = "mkdir " + to_string(p) + "_m";
    system(dir.c_str());

    for (mp = 1; mp <= 6; mp++)
    {

        fstream file1;
        ll n = power(p, mp);
        string name1 = to_string(p) + "_m/" + to_string(n) + "_cr_analysis.html";
        string name2 = to_string(p) + "_m/" + to_string(n) + "_qr_analysis.html";
        string name3 = to_string(p) + "_m/" + to_string(n) + "_6r_analysis.html";

        string n1 = "../" + name1;
        string n2 = "../" + name2;
        string n3 = "../" + name3;

        file1.open(to_string(p) + "_m/" + to_string(n) + "_cr_qr_analysis.html", ios::out);

        file1 << "<!DOCTYPE html>\n";
        file1 << "<html>\n";
        file1 << "<head>\n";

        file1 << "<title>Analysis of Cubic Residue and Quadratic Residue over Z<sub>" + to_string(n) + "</sub></title>\n";

        file1 << "<style>\n";
        file1 << "body { font-family: Arial, sans-serif; justify-content: center; align-items: center; }\n";
        file1 << ".container { display: flex; flex-direction: column; align-items: center; }\n";
        file1 << ".container button { background-color: #dae65a; color: white; padding: 10px 20px; border: none; cursor: pointer; margin: 5px; width: 310px; }\n";
        file1 << ".container button:hover { background-color: #dceb3d; }\n";
        file1 << "a { text-decoration: none; color: #333;} \n a:hover {color: #007bff;}";
        file1 << "</style>\n";
        file1 << "</head>\n";
        file1 << "<body>\n";
        file1 << "<center>\n";
        file1 << "<h4>N = " << n << "</h4>\n";
        file1 << "</center>\n";
        file1 << "<div class=\"container\">\n";
        file1 << "<button><a href=\"" + n2 + "\">Quadratic Residue analysis</a></button>\n";
        file1 << "<button><a href=\"" + n1 + "\">Cubic Residue analysis</a></button>\n";
        file1 << "<button><a href=\"" + n3 + "\">sixth Residue analysis</a></button>\n";
        file1 << "</div>\n";
        file1 << "</body>\n";
        file1 << "</html>\n";
        file1.close();

        fstream file2, file3, file4;

        file2.open(name1, ios::out);
        file3.open(name2, ios::out);
        file4.open(name3, ios::out);

        ///////////////////////////////////////////CR analysis///////////////////////////////////////////////
        file2 << "<!DOCTYPE html>\n";
        file2 << "<html>\n";
        file2 << "<head>\n";
        file2 << "<title>Analysis of Cubic Residue over Z<sub>" + to_string(n) + "</sub></title>\n";
        file2 << "<style>\n";
        file2 << "body { font-family: Arial, sans-serif; align-items: center; justify-content: flex-start; height: 100vh; margin: 0; }\n";
        file2 << "table { width: 50%; margin: 0 auto; border-collapse: collapse; }\n";
        file2 << "th, td { padding: 10px; border: 1px solid #333; }\n";
        file2 << "th { background-color: #faff75; }\n";
        file2 << "tr:nth-child(even) { background-color: #f2f2f2; }\n";
        file2 << "tr:hover { background-color: #ddd; }\n";
        file2 << ".container { text-align: center; }\n";
        file2 << "</style>\n";
        file2 << "</head>\n";
        file2 << "<body>\n";
        file2 << "<div class=\"container\">\n";

        file2 << "<h3>Analysis of Cubic Residue over Z<sub>" + to_string(n) + "</sub></h3>\n";
        file2 << "<table>\n";

        file2 << "<tr><th>Cubic residue (CR) </th><th>x such that x<sup>3</sup> = CR mod " + to_string(n) + "</th><th>Occurance</th></tr>\n";

        // map to store cubic residue with its occurance and x such that x^3 = CR mod p^m
        map<ll, pair<ll, vector<ll>>> m_cr;
        for (ll i = 0; i < n; i++)
        {
            ll x = (i % n * (i % n * i % n) % n) % n;
            m_cr[x].first++;
            m_cr[x].second.push_back(i);
        }

        for (auto i = m_cr.begin(); i != m_cr.end(); i++)
        {
            file2 << "<tr><td>" << i->first << "</td><td>";
            
            for (auto j = i->second.second.begin(); j != i->second.second.end() - 1; j++)
                file2 << *j << ", ";
            
            file2 << *(i->second.second.end() - 1) << "</td><td>" << i->second.first << "</td></tr>\n";
        }
        file2 << "</table>\n";
        file2 << "</div>\n";

        // map to store occurance of cubic residue with cubic residue
        map<ll, vector<ll>> m_occ_cr;
        for (auto i = m_cr.begin(); i != m_cr.end(); i++)
            m_occ_cr[i->second.first].push_back(i->first);
        file2 << "<div class=\"container\" style=\"margin-top: 1%;\">\n";
        file2 << "<table>\n";
        file2 << "<tr><th>Occurance</th><th>Cubic residues</th><th>Total CRs with same occurance</th></tr>\n";

        for (auto i = m_occ_cr.begin(); i != m_occ_cr.end(); i++)
        {
            file2 << "<tr><td>" << i->first << "</td><td>";

            for (ll j = 0; j < i->second.size() - 1; j++)
                file2 << i->second[j] << ", ";
            file2 << i->second[i->second.size() - 1] << "</td><td>" << i->second.size() << "</td></tr>\n";
        }
        file2 << "</table>\n";
        file2 << "</div>\n";
        ll cr_cnt = m_cr.size();
        
        file2 << "<center><h3>Total Cubic Residues = " << cr_cnt << "</h3></center>\n";
        file2 << "</body>\n";
        file2 << "</html>\n";
        file2.close();
        ///////////////////////////////////////////QR analysis///////////////////////////////////////////////
        map<ll, pair<ll, vector<ll>>> m_qr;
        for (ll i = 0; i < n; i++)
        {
            ll x = ((i % n) * (i % n)) % n;
            m_qr[x].first++;
            m_qr[x].second.push_back(i);
        }
        ll qr_cnt = m_qr.size();

        map<ll, vector<ll>> m_occ_qr;
        for (auto i = m_qr.begin(); i != m_qr.end(); i++)
            m_occ_qr[i->second.first].push_back(i->first);
        
        file3 << "<!DOCTYPE html>\n";
        file3 << "<html>\n";
        file3 << "<head>\n";

        file3 << "<title>Analysis of Quadratic Residue over Z<sub>" + to_string(n) + "</sub></title>\n";
        file3 << "<style>\n";
        file3 << "body { font-family: Arial, sans-serif; align-items: center; justify-content: flex-start; height: 100vh; margin: 0; }\n";
        file3 << "table { width: 50%; margin: 0 auto; border-collapse: collapse; }\n";
        file3 << "th, td { padding: 10px; border: 1px solid #333; }\n";
        file3 << "th { background-color: #faff75; }\n";
        file3 << "tr:nth-child(even) { background-color: #f2f2f2; }\n";
        file3 << "tr:hover { background-color: #ddd; }\n";
        file3 << ".container { text-align: center; }\n";
        file3 << "</style>\n";
        file3 << "</head>\n";
        file3 << "<body>\n";
        file3 << "<div class=\"container\">\n";

        file3 << "<h3>Analysis of Quadratic Residue over Z<sub>" + to_string(n) + "</sub></h3>\n";
        file3 << "<table>\n";

        file3 << "<tr><th>Quadratic residue (QR) </th><th>x such that x<sup>2</sup> = QR mod " + to_string(n) + "</th><th>Occurance</th></tr>\n";

        for (auto i = m_qr.begin(); i != m_qr.end(); i++)
        {
            file3 << "<tr><td>" << i->first << "</td><td>";
           
            for (auto j = i->second.second.begin(); j != i->second.second.end() - 1; j++)
                file3 << *j << ", ";
            
            file3 << *(i->second.second.end() - 1) << "</td><td>" << i->second.first << "</td></tr>\n";
        }
        file3 << "</table>\n";
        file3 << "</div>\n";
        file3 << "<div class=\"container\" style=\"margin-top: 1%;\">\n";
        file3 << "<table>\n";
        file3 << "<tr><th>Occurance</th><th>Quadratic residues</th><th>Total QRs with same occurance</th></tr>\n";
        for (auto i = m_occ_qr.begin(); i != m_occ_qr.end(); i++)
        {
            file3 << "<tr><td>" << i->first << "</td><td>";
            
            for (auto j = i->second.begin(); j != i->second.end() - 1; j++)
                file3 << *j << ", ";
            
            file3 << *(i->second.end() - 1) << "</td><td>" << i->second.size() << "</td></tr>\n";
        }
        file3 << "</table>\n";
        file3 << "</div>\n";
        file3 << "<center><h3>Total Quadratic Residues = " << qr_cnt << "</h3></center>\n";
        file3 << "</body>\n";
        file3 << "</html>\n";
        file3.close();
        ///////////////////////////////////////////6R analysis///////////////////////////////////////////////
        map<ll, pair<ll, vector<ll>>> m_6r;
        for (ll i = 0; i < n; i++)
        {
            ll x = mod_pow(i, 6, n);
            m_6r[x].first++;
            m_6r[x].second.push_back(i);
        }
        ll sixr_cnt = m_6r.size();
        map<ll, vector<ll>> m_occ_6r;
        for (auto i = m_6r.begin(); i != m_6r.end(); i++)
            m_occ_6r[i->second.first].push_back(i->first);

        
        file4 << "<!DOCTYPE html>\n";
        file4 << "<html>\n";
        file4 << "<head>\n";

        file4 << "<title>Analysis of 6th Residue over Z<sub>" + to_string(n) + "</sub></title>\n";
        file4 << "<style>\n";
        file4 << "body { font-family: Arial, sans-serif; align-items: center; justify-content: flex-start; height: 100vh; margin: 0; }\n";
        file4 << "table { width: 50%; margin: 0 auto; border-collapse: collapse; }\n";
        file4 << "th, td { padding: 10px; border: 1px solid #333; }\n";
        file4 << "th { background-color: #faff75; }\n";
        file4 << "tr:nth-child(even) { background-color: #f2f2f2; }\n";
        file4 << "tr:hover { background-color: #ddd; }\n";
        file4 << ".container { text-align: center; }\n";
        file4 << "</style>\n";
        file4 << "</head>\n";
        file4 << "<body>\n";
        file4 << "<div class=\"container\">\n";

        file4 << "<h3>Analysis of sixth Residue over Z<sub>" + to_string(n) + "</sub></h3>\n";
        file4 << "<table>\n";
        file4 << "<tr><th>sixth residue (6R) </th><th>x such that x<sup>6</sup> = 6R mod " + to_string(n) + "</th><th>Occurance</th></tr>\n";
        for (auto i = m_6r.begin(); i != m_6r.end(); i++)
        {
            file4 << "<tr><td>" << i->first << "</td><td>";
            
            for (auto j = i->second.second.begin(); j != i->second.second.end() - 1; j++)
                file4 << *j << ", ";
            
            file4 << *(i->second.second.end() - 1) << "</td><td>" << i->second.first << "</td></tr>\n";
        }
        file4 << "</table>\n";
        file4 << "</div>\n";
        file4 << "<div class=\"container\" style=\"margin-top: 1%;\">\n";
        file4 << "<table>\n";
        file4 << "<tr><th>Occurance</th><th>sixth residues</th><th>Total 6Rs with same occurance</th></tr>\n";
        for (auto i = m_occ_6r.begin(); i != m_occ_6r.end(); i++)
        {
            file4 << "<tr><td>" << i->first << "</td><td>";
            
            for (auto j = i->second.begin(); j != i->second.end() - 1; j++)
                file4 << *j << ", ";
            
            file4 << *(i->second.end() - 1) << "</td><td>" << i->second.size() << "</td></tr>\n";
        }
        file4 << "</table>\n";
        file4 << "</div>\n";
        file4 << "<center><h3>Total sixth Residues = " << sixr_cnt << "</h3></center>\n";
        file4 << "</body>\n";
        file4 << "</html>\n";
        file4.close();
    }
    return 0;
}