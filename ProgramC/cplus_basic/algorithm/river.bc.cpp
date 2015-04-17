#include <stdio.h>
#include <iostream>
using namespace std;


bool PsSet[4][4], turn;//PsSet是允许状态集 turn表示由此岸到彼岸或由彼岸到此岸
int PdSet[5][2]={{0,1},{0,2},{2,0},{1,0},{1,1}};//PdSet表示允许决策集
int Record[12][2],pos=1;//Record记录已走过的状态 pos当前位置
bool Check[12];//Check记录走过状态的turn值
int num=0;

void InitPsSet(){  //初始化PsSet
   for(int i=0;i<4;i++)
       for(int j=0;j<4;j++)
           if(i==0||i==3||i==j);
           else PsSet[i][j]=true;
}

int IsRepeat(int px,int py,bool turn) //检验该状态是否被用过
{ 
    for(int i=0;i<pos;i++)
        if(px==Record[i][0]&&py==Record[i][1]&&turn==Check[i])
              return true;
	else
              return false;
}

void CrossRiver(int px,int py)  //过河子函数
{ 
    int i;
    if(px==0&&py==0)
    { 
       num++;
       printf("第 %d 种方案",num);
       printf("\n");
       for(i=0;i<pos;i++)
       {  
           printf("(%d,%d) 到 (%d,%d)",Record[i][0],Record[i][1],3-Record[i][0],3-Record[i][1]);
           printf("\n");
       }
     }else if(turn==false)
     { 
        for(i=0;i<5;i++)
        {
	    px=px-PdSet[i][0];
            py=py-PdSet[i][1];
            if(px<0||py<0||px>3||py>3||PsSet[px][py]||IsRepeat(px,py,turn))
            { 
	         px=px+PdSet[i][0];
                 py=py+PdSet[i][1];
                 continue;
            }
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
        }
      }else {
            for(i=0;i<5;i++)
            {
	        px=px+PdSet[i][0];
                py=py+PdSet[i][1];
                if(px<0||py<0||px>3||py>3||PsSet[px][py]||IsRepeat(px,py,turn))
                { 
		    px=px-PdSet[i][0];
                    py=py-PdSet[i][1];
                    continue;
                }
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
            }
       }
}

int main(){  //主函数
int px=3,py=3;
InitPsSet();
Check[0]=true;
Record[0][0]=3;
Record[0][1]=3;
cout<<px<<endl;
//CrossRiver(px,py);
}

