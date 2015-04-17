#include <stdio.h>
#include <iostream>
using namespace std;


bool PsSet[4][4], turn;// PsSet: possible state set, turn: forward or backward
int PdSet[5][2]={{0,1},{0,2},{2,0},{1,0},{1,1}};// PdSet: possible disicion set
int Record[12][2],pos=1;//Record: record the road, pos: current position
bool Check[12];//Check 
int num=0;

void InitPsSet(){  //Initial PsSet, only valid states were set to true
   cout<<"People"<<"\t"<<"Wolf"<<"\t"<<"live or die"<<endl;
   for(int i=0;i<4;i++)
       for(int j=0;j<4;j++)
           if(i==0||i==3||i==j){
              PsSet[i][j]=true;
              cout<<i<<"\t"<<j<<"\t"<<"true"<<endl;
           }
}

int IsRepeat(int px,int py,bool turn) //检验该状态是否被用过
{ 
    for(int i=0;i<pos;i++)
        if(px==Record[i][0]&&py==Record[i][1]&&turn==Check[i])
              return true;
	else
              return false;
}

void CrossRiver(int px,int py)  // cross river with boat, px: number of people, py: number of wolf
{ 
    int i;
    if(px==0&&py==0) // mission complete
    { 
       num++;
       //printf("strategy num: %d",num);
       //printf("\n");
       for(i=0;i<pos;i++)
       {  
           //printf("(%d,%d) to (%d,%d)",Record[i][0],Record[i][1],3-Record[i][0],3-Record[i][1]);
           //printf("\n");
           cout<<"test"<<endl;
       }
    }else if(turn==false) // toward the target side
    { 
        for(i=0;i<5;i++)
        {
	    px=px-PdSet[i][0];
            py=py-PdSet[i][1];
            cout<<"Forw: "<<px<<"\t"<<py<<endl;
            if(px<0||py<0||px>3||py>3||PsSet[px][py]||IsRepeat(px,py,turn))
            { 
	         px=px+PdSet[i][0];
                 py=py+PdSet[i][1];
                 continue;
            }
            cout<<"Forw2: "<<px<<"\t"<<py<<endl;
            Record[pos][0]=px;
            Record[pos][1]=py;
            Check[pos]=turn;
            pos++;
            turn=!turn;
            CrossRiver(px,py);
            pos--;
            turn=!turn;
            px=px+PdSet[i][0];
            py=py+PdSet[i][1];
            cout<<"Forw3: "<<px<<"\t"<<py<<endl;
        }
    }else { // backward to the source side
            for(i=0;i<5;i++) // foreach PdSet strategy
            {
	        px=px+PdSet[i][0];
                py=py+PdSet[i][1];
                cout<<"Back: "<<px<<"\t"<<py<<endl;
                if(px<0||py<0||px>3||py>3||PsSet[px][py]||IsRepeat(px,py,turn)) //
                { 
		    px=px-PdSet[i][0];
                    py=py-PdSet[i][1];
                    continue;
                }
                cout<<"Back2: "<<px<<"\t"<<py<<endl;
                Record[pos][0]=px;
                Record[pos][1]=py;
                Check[pos]=turn;
                pos++;
                turn=!turn;
                CrossRiver(px,py);
                pos--;
                turn=!turn;
                px=px-PdSet[i][0];
                py=py-PdSet[i][1];
                cout<<"Back3: "<<px<<"\t"<<py<<endl;
            }
    }
}

int main(){  //主函数
int px=3,py=3;
InitPsSet();
Check[0]=true;
Record[0][0]=3;
Record[0][1]=3;
//cout<<px<<endl;
CrossRiver(px,py);
}

