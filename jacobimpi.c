#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <mpi.h>
#include <string.h>
int size, rank;

double* linspace(double a, double b, int n){
    //this subroutine mimicks the matlab function linspace
    double c; int i;
    double *u = malloc(n*sizeof(double));
    c = (b - a)/(n - 1);
    for(i=0; i<n-1; ++i)
        u[i] = a + i*c;
    u[n-1] = b;
    return u;
}

double func(double x, double y){
    //forcing function
    return exp(x+y)*((x*x+3*x)*(y*y-y)+(y*y+3*y)*(x*x-x));
}

double sol(double x, double y){
    //Exact Solution (This won't equal the fully converged solution as there is a 2nd order error)
    return exp(x+y)*(x*x-x)*(y*y-y);
}

double finderr(int m, double a[m][m], double b[m][m]){
    int i,j;
    double maxi = -1e13;
    for(i=1; i<m-1; i++){
        for(j=1; j<m-1; j++)
            maxi = (fabs(a[i][j]-b[i][j])>maxi)?fabs(a[i][j]-b[i][j]):maxi;}
    return maxi;
}

int main(int argc, char **argv){
    double threshold; int jet0=0, jet1=0, iterCount;
    if (!(argc==3)){
        if (rank==0) printf("Error Wrong Argument Count. Assuming threshold 1e-4 and iteration limit 4000.\n");
        threshold = 1e-4; iterCount = 4000;
    }
    else{
      jet0 = sscanf(argv[1], "%i" , &iterCount);
      jet1 = sscanf(argv[2], "%lg" , &threshold);
    }
    if ((jet0==-1) || (jet1==-1)){
        if (rank==0) printf("Error Wrong Argument Type\n"); return -1;
    }
    int i,j,p,q,k=0,xlim,ylim;
    MPI_Init(&argc, &argv); //Open global comm network
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    p = rank%4;
    q = rank/4;
    
    double *x = linspace(0+p/4.0,0.25+1.0/256+p/4.0,66);
    double *y = linspace(0+q/4.0,0.25+1.0/256+q/4.0,66);
    
    double u[66][66], f[66][66], v[66][66], s[66][66], b[66][66], west[66], westrecv[66], east[66], eastrecv[66],
            north[66], northrecv[66], south[66], southrecv[66];
    for (j = 0; j < 66; j++){
        for (i = 0; i < 66; i++){
            u[i][j] = 0.0 * rank;
            v[i][j] = 0.0;
            f[i][j] = func(x[i],y[j]);
            s[i][j] = sol (x[i],y[j]);
            b[i][j] = 0.0;
        }
    }
    double error, gerror, diff, gdiff;
    while(1){
      if((k%10)==0){
          if (k!=0){
            /*
             ######      ###    ######## ##     ## ######## ########        ######## ########  ########      
            ##    ##    ## ##      ##    ##     ## ##       ##     ##       ##       ##     ## ##     ##     
            ##         ##   ##     ##    ##     ## ##       ##     ##       ##       ##     ## ##     ##     
            ##   #### ##     ##    ##    ######### ######   ########        ######   ########  ########      
            ##    ##  #########    ##    ##     ## ##       ##   ##         ##       ##   ##   ##   ##       
            ##    ##  ##     ##    ##    ##     ## ##       ##    ##        ##       ##    ##  ##    ##  ## 
             ######   ##     ##    ##    ##     ## ######## ##     ##       ######## ##     ## ##     ## ors 
            */
            diff = finderr(66,u,s);
            MPI_Allreduce(&diff, &gdiff, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
            error = finderr(66,b,u);
            MPI_Allreduce(&error, &gerror, 1, MPI_DOUBLE, MPI_MAX,  MPI_COMM_WORLD);
            for(i=1; i<65; i++){
              for(j=1; j<65; j++)
                b[i][j] = u[i][j];
            }
            
            MPI_Barrier(MPI_COMM_WORLD);
            //if (rank==0) printf("%i\t%e\t%e\n",k,gerror,gdiff);
            if ((k >= iterCount) && (gerror<threshold)){
              if (rank==0) printf("d = %e\tc = %e\n\n",gerror,gdiff);
              break;
            }
          }
      }
      k++; 
        /*
        ########     ###    ########    ###          ######## ##     ##  ######  ##     ## 
        ##     ##   ## ##      ##      ## ##         ##        ##   ##  ##    ## ##     ##       
        ##     ##  ##   ##     ##     ##   ##        ##         ## ##   ##       ##     ##       
        ##     ## ##     ##    ##    ##     ##       ######      ###    ##       #########   
        ##     ## #########    ##    #########       ##         ## ##   ##       ##     ##       
        ##     ## ##     ##    ##    ##     ##       ##        ##   ##  ##    ## ##     ## ###       
        ########  ##     ##    ##    ##     ##       ######## ##     ##  ######  ##     ##ange
        */ 
        MPI_Barrier(MPI_COMM_WORLD);
        for (j = 0; j < 66; j++){
            west [j] = u[64][j];
            east [j] = u[1] [j];
            north[j] = u[j][64];
            south[j] = u[j] [1];
        }
        //MPI_Barrier(MPI_COMM_WORLD);
        
        //Send West Data to east of west grid
        if ((p==0) || (p==1) || (p==2))
            MPI_Send(west,     66, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD);
        if ((p==1) || (p==2) || (p==3)){
            MPI_Recv(westrecv, 66, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (j = 0; j < 66; j++)
                u[0][j] = westrecv[j];
        }
        
        //MPI_Barrier(MPI_COMM_WORLD);
        
        //Send East Data of grid to west of eastern grid
        if ((p==1) || (p==2) || (p==3))
            MPI_Send(east,     66, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD);
        if ((p==0) || (p==1) || (p==2)){
            MPI_Recv(eastrecv, 66, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (j = 0; j < 66; j++)
                u[65][j] = eastrecv[j];        
        }
        
        //MPI_Barrier(MPI_COMM_WORLD);
        
        //Send north Data to south of northern grid
        if ((q==0) || (q==1) || (q==2))
            MPI_Send(north,     66, MPI_DOUBLE, rank+4, 0, MPI_COMM_WORLD);
        if ((q==1) || (q==2) || (q==3)){
            MPI_Recv(northrecv, 66, MPI_DOUBLE, rank-4, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (i = 0; i < 66; i++)
                u[i][0] = northrecv[i];
        }
        
        //MPI_Barrier(MPI_COMM_WORLD);
        
        //Send south Data of grid to north of southern grid
        if ((q==1) || (q==2) || (q==3))
            MPI_Send(south,     66, MPI_DOUBLE, rank-4, 0, MPI_COMM_WORLD);
        if ((q==0) || (q==1) || (q==2)){
            MPI_Recv(southrecv, 66, MPI_DOUBLE, rank+4, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);    
            for (i = 0; i < 66; i++)
                u[i][65] = southrecv[i];
        }
        
        MPI_Barrier(MPI_COMM_WORLD); // syncronize
        /*
              ##    ###     ######   #######  ########  ####       #### ######## ######## ########  
              ##   ## ##   ##    ## ##     ## ##     ##  ##         ##     ##    ##       ##     ## 
              ##  ##   ##  ##       ##     ## ##     ##  ##         ##     ##    ##       ##     ## 
              ## ##     ## ##       ##     ## ########   ##         ##     ##    ######   ########  
        ##    ## ######### ##       ##     ## ##     ##  ##         ##     ##    ##       ##   ##   
        ##    ## ##     ## ##    ## ##     ## ##     ##  ##         ##     ##    ##       ##    ##   ###
         ######  ##     ##  ######   #######  ########  ####       ####    ##    ######## ##     ##ation
        */
        xlim = (p==3)?64:65;
        ylim = (q==3)?64:65;
        for(i=1; i<xlim; i++){
            for(j=1; j<ylim; j++){
                v[i][j] = 0.25*(u[i+1][j] + u[i-1][j] + u[i][j+1] + u[i][j-1] - f[i][j]/65536.0);
            }
        }
        

        // Backstore v on u
        for(i=1; i<xlim; i++){
            for(j=1; j<ylim; j++){
                u[i][j] = v[i][j];
            }
        }
        
    }
    
    MPI_Finalize(); //Close comm channel
    return 0;
}
/* Coder: Adil Ansari
ssso++++++ooooyo+ooo++++/++////oososo/:::::::::--::--------------------::::::-::--:::::
ooooooo++ooo++++++o+++///+++++/+ossyyso+//:://::::::::::----------------:::::::::::::::
oooooo+++oo++++//:::/:/++////+/++++ydhy+ooosssssyyssyhdy+/:-------------:::::::::::::::
ooo++o+++oo++++++//////+o+//++++///ohddssyhdmmNNNmdddddmmdysss+/::------:::::::::::::::
oo++//++oo+oooo++++/////+//++/+o+//ydNNNmmNNNmmmmNmmmmmmNNNNNNmh+:::::::::::::::::::://
+++///////////++osso++/////++++oooyNMMMMNNNNNmmmNNNmmmmNMMMMMMNNmo/::::::::::::::::////
+ooo+++++o////+ossyyyo+/+/////+osymMMMMMMNNNNNMMMMNNNNNNNNMMMMMMMmyo:::::::::::::::////
sssoooo++o///++oso+shysoooo+++osydNMMMMMMMMMMMMMMNmmmmmNmmmNMMMMMMNmh+::::::::::///////
soooooo+////++/:/+++++osysyyhhhhdmNMMMMNMMMMMMMMMNNNNNNNNNNmMMMMMMMMNd/::/:::://///////
++++++++///////++++//+++++oydNNmmmNNmdddmmmmmmNmmdNNMMMMMMNNMMMMMMMMMNs////////////////
+++++/++++++/://++/////////ymMNNmddhs+//oooooo++++ydmNMMMMNNNNMMMMMMMMNs+////////////++
+++/::///////////////////+ohNMNmmmmh+:-.---..``.-:+yhdmNNNNNNNMMMMMMMMMho//////////++++
o+///////////:/++++++//::-+hNNNNNNNy/-.````    ```.:+oyhdmmNNNNMMMMMMMds/:://////++++++
+///:---:/:::/++////o+++/:odNNMNNmd/:..`````````````..-:+oyhdmmNNMMMMMo/::::///++++++++
++++/::::/:-:://///+oosso/sdNNMNho/----....`````````````.-:+yhddmMMMMN+//+/+/:://++oo++
//////+////////:::/+//::::ohNMNy+:-://+osssoo/:--...-...---:+syhdMMMMm///+ss+/---:::///
//////ooooo++/:-:::://:::/+ymMhyyso+++oydmmmmmdyo--+ossyyyssosyyhNMMNy:://o++/::::--:::
///+/+sssoo+/:::::/+++///ssyhh+hdh/oyhhdmmmmmmmmmyddddmNNNmmdmmdhNNNdo::://++/::/::::--
/++++ooooooo++///++++++///:+s+-+ss--/+ossyyyyyhyo/hdddmNNNNmmdddNmmhso/////+////////:--
+ooo++++///+oo/+//+++++//..:+:.-+o:...-:::---+/-``/ssossyyhdmdhhmdmdy+++++++/////:://--
ooooooo+/////++//+osoo++/...::...-:......--::.`` `./+/--::/+oooo+hmdysssoos+/://::://::
//++ooo+//://///://++++//.``.:-........-/+:-```  `.-/+/:--:://::/hdhsosooos+//+o////+//
//++//+++++++//o+////+///-.`.::--.....-/o:--.....-::/+o/::::://+odyo//:::////++o//+/://
///+//oo+/////://////++++/:-:///::::-:/+--/sysoosyhhhs+::://+oosss+/::::::://+++::///++
:///+so/:://///+////////://////////:://:-.-:/yhddmmddhs///++sssss+///:-:/:-:::::-:/+s+/
///++///+//////+++++//:::::::++/////:/+syo+///osyyhhddy/+osyyys+::::://+//::::::::/+y+/
/+++///+++////+ooooo++////:::+++///////://::::/+oydmNmdoosyyyyo/:/++++++++/:::---:++o/:
oo++/////:///++soooso+++/////+oo+++//::--------:/ohhdhysyyhhyo+//:/ossshdhy//::----::-:
ooo++/::::////++///++/////////ssso+/::---::+oyyyysoooooshhhhs++/:-:/osshdho/:---:---:-:
////////:://////////////:::::-oyyyo+/:--..---::::/::/+ohhhs+ohhy/:/++/+++/::::::::-----
::::-:::::-:///+++/+ooo++/:--.-oyhhys+/--.....---::+sydds+//+oss+/////:///:::://::::::-
::://:::::::::/+++++:--..---...-+ydddhyo++//+++oooyhdddy+/:////////++//+++/:///+o++////
`.-://///++++++ooo+-`````.......-+yddddddddddmmmmmmmmddhhyo:::://:::::/+o+/://///:::::/
:::/+oooosooooosss+.``````.......-:shhdddmmmmmmmmmdddddmNNmy//////////:::::////::::::::
sssssssssoooooossso:`````.........-/oyhhdddmmmmmdddddddmNNNNmdyo///////::::+++////:::::
ssssssssssoooossssso.```````......--/oydddmmmdddddddddmNNMMMMMNNdo+////:::::////::--:::
ssssssssssossssoossso:.``````.....-:/+shddddddddddddddNMMMMMMNNmdhso++/:-------::-----:
ssssssssssosssssoosssso/-````......--:/+hdddddddddmmNMMMMMMNmhysooooo++/:---.......-..-
sssssssssssssssssssssssss+/:--...----.../shdmmmmNNMMMMNMNNmhsoooooooooo++:---....----..
ssssssssssssssssssssssssyyyysso++//:::::/syhmNNMMMMMMNNmdysooooooooooo+++++:----:///:--
*/

