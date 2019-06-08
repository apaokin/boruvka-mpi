#include <iostream>
#include <map>
#include <omp.h>

using namespace std;

int main()
{
  #pragma omp parallel for
  for(int k=0; k < 5;k++){
    printf()
  }
  map <string,int> myMap = {{ "Beta", 2 }, ///явная инициализация map
                            { "Alpha", 1 },
                            { "Gamma", 3 }};

  ///присвоение элементам map новых значений
   myMap.at("Beta") = 0;
   myMap.at("Alpha") = 233;
   myMap.at("Gamma") = -45;
   // cout << myMap.at("FFF")<<endl;
   cout<<( myMap.find( "FFF" ) != myMap.end())<<endl;
  cout << "myMap contains:\n";
  for(auto it = myMap.begin(); it != myMap.end(); ++it)
  {
      cout << it->first << " : " << it->second << endl;///вывод на экран
  }

  multimap <char,int> myMultimap;///объявили multimap

   ///заполняем myMultimap
   myMultimap.insert ( pair<char,int>('q',111) );
   myMultimap.insert ( pair<char,int>('u',201) );
   myMultimap.insert ( pair<char,int>('h',301) );

   cout << "\nmyMultimap contains:\n";
   for (auto it = myMultimap.begin(); it != myMultimap.end(); ++it)
   {
      cout << it->first << " : " << it->second << endl;
   }

   myMap.clear();
   myMultimap.clear();

   ///новая инициализация myMap
   myMap = {{ "Mike", 40 },
            { "Walle", 999 },
            { "Cloude", 17 }};

   ///новая инициализация myMultimap
   myMultimap.insert ( pair<char,int>('q',222) );
   myMultimap.insert ( pair<char,int>('u',223) );
   myMultimap.insert ( pair<char,int>('h',221) );

   auto itMap = myMap.begin();///создаем итератор на начало myМap
   auto itMultimap = myMultimap.begin();///создаем итератор на начало myMultimap
   cout << "\nmyMap after clear contains: \t myMultimap after clear contains:\n";

   ///вывод на экран myMap и myMultimap
   for(itMap = myMap.begin(),itMultimap = myMultimap.begin(); itMultimap != myMultimap.end(); itMap++,itMultimap++)
   {
       cout << "\t" <<itMap->first << " : " << itMap->second << "\t\t\t\t" << itMultimap->first << " : " << itMultimap->second << endl;
   }
   return 0;
}
