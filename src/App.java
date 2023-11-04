import Jama.*;
import javax.swing.JFrame;
//import org.math.plot.Plot2DPanel;
import org.math.plot.Plot3DPanel;
public class App {
    static int n=3;
    static int m=15;
    static void MatToVec(double[][] mat,double[] V1 , double[] V2, double[] V3){
        
        for(int i=0;i<m;i++){
            V1[i]=mat[i][0];
            V2[i]=mat[i][1];
            V3[i]=mat[i][2];
            
        }
    }
    static void g3D(double[] V1 , double[] V2, double[] V3){
        Plot3DPanel plot = new Plot3DPanel(); 
        plot.addScatterPlot("Indv3D", V1, V2,V3);
        JFrame frame = new JFrame("Les individus");
        frame.setContentPane(plot);
        frame.setSize(850,550);
        frame.setVisible(true);
        
    
    }
    static void affichermat(double tab[][]){
        System.out.println("\n\n\t************************* MATRICE Initial *************************\n");
        for(int i = 0; i <m; i++) {
			for(int j = 0; j < n; j++) {
				System.out.printf("\t%.1f\t", tab[i][j]);
			}
			System.out.println();
		}
    }

    static void centrée(double[][] tab,double[][] matC){
        System.out.println("\n\n\t************************* MATRICE CENTREE *************************\n");
        double moy=0;
        for(int j = 0; j < n; j++) {
			for(int i = 0; i <m; i++) {
				moy = moy + tab[i][j];
			}
			moy = moy / m;
			for(int i = 0; i <m; i++) {
				matC[i][j] = tab[i][j] - moy;
			}
			moy = 0;
		}
        System.out.println("\t x1\t\t\t x2\t\t\t x3\n");
        System.out.println("\t______________________________________________________\n\n");
		for(int i = 0; i <m; i++) {
			for(int j = 0; j < n; j++) {
				System.out.printf("\t%f\t", matC[i][j]);
			}
			System.out.println("\n");
		}
        
    }

    static void reduit(double[][] matC,double[][] matR){
        System.out.println("\n\n\t*************************MATRICE REDUITE*************************\n");
        double[] sigma = new double[n];
        for(int j = 0; j < n; j++) {
			sigma[j] = 0;
			for(int i = 0; i <m; i++) {
				sigma[j] = sigma[j] + matC[i][j] * matC[i][j];
			}
			sigma[j] = sigma[j] / m;
			for(int i = 0; i <m; i++) {
				matR[i][j] = (double)(matC[i][j] / Math.sqrt(sigma[j]));
			}
		}
        System.out.println("\t x1\t\t\t x2\t\t\t x3\n");
        System.out.println("\t______________________________________________________\n\n");
		for(int i = 0; i <m; i++) {
			for(int j = 0; j < n; j++) {
				System.out.printf("\t%f\t", matR[i][j]);
			}
			System.out.println("\n");
		}
    }
    static void Normees(double[][] matR,double[][] matN){
        System.out.println("\n\n\t************************* MATRICE NORMEES *************************\n");
        for(int i=0;i<m;i++){
            for(int j=0;j<n;j++){
                matN[i][j]= (double) (matR[i][j] / Math.sqrt(m));
            }
        }
        System.out.println("\t x1\t\t\t x2\t\t\t x3\n");
        System.out.println("\t______________________________________________________\n\n");
        for(int i=0;i<m;i++){
            for(int j=0;j<n;j++){
                System.out.printf("\t%f\t", matN[i][j]);
            }
            System.out.println("\n");
        }
    }
    static void transpose(double[][] matN, double [][] matNT){
        System.out.println("\n\n\t************************* MATRICE NORMEES Transposé *************************\n");
        for(int i=0;i<m;i++){
            for(int j=0;j<n;j++){
                matNT[j][i]=matN[i][j];
            }
        }
        for(int i=0;i<3;i++){
            for(int j=0;j<m;j++){
                System.out.printf("\t%f\t", matNT[i][j]);
            }
            System.out.println("\n______________________________________________________________________________________________________________________________________________________________\n");
            System.out.println("\n");
        }
    }
    static void CORRELATION(double[][] matNT,double[][] matN,double[][] corr){
        System.out.println("\n\n\t************************* MATRICE CORRELATION *************************\n\n");
        for (int i=0;i<n;i++){
            for (int j=0;j<n;j++){
                    double s=0;
                for (int k=0;k<m;k++){
                    s+=matNT[i][k]*matN[k][j];
                }
                corr[i][j]=s;
            }
        }
        for(int i=0;i<n;i++){
            for(int j=0;j<n;j++){
                System.out.printf("\t%f\t",corr[i][j]);
            }
            System.out.println("\n");
        }
    }
    static void vecteurPr(double[][] corr,double[][] vep2){
        System.out.println("\n\n\t *************************Les vecteurs Propre1 *************************\n\n");
        Matrix m1= new Matrix(corr);
        double[][] vep3=new double[n][n];
        vep3= m1.eig().getV().getArray();
        for(int i=0 ; i<n ;i++){
            for(int j=0;j<n ;j++){
                vep2[i][j]=vep3[i][j];
            System.out.print(vep2[i][j]);
            System.out.print("\t\t");
            }
            System.out.println("\n");
        }
    }
    static void vecteurPr2(double[][] vep1,double[][] vep2){
        System.out.println("\n\n\t *************************Les vecteurs Propre2 *************************\n\n");
        for(int i=0;i<n ;i++){
                for(int j=0;j<n ;j++){
                    vep1[j][i]=vep2[i][j];
                }
        }
            double p1;
            int  x=0,m=2;
            for(int j=0;j<n;j++){
                p1=vep1[x][j];
                vep1[x][j]=vep1[m][j];
                vep1[m][j]=p1;
            }
    
        for(int i=0 ; i<n ;i++){
            for(int j=0;j<n ;j++){
                System.out.print(vep1[i][j]);
                System.out.print("\t\t");
            }
            System.out.println("\n");
        }
    }
    static void vecteurPrTr(double [][] vep1){
        System.out.println("\n\n\t************************* Les vecteurs tri Propre ************************* \n\n");
        double p;
        for(int i=0;i<n;i++){
            for(int j=0;j<n ; j++){
                for(int k=0; k<n ; k++){
                    if(vep1[i][k]<vep1[i][j]){
                        p=vep1[i][j];
                        vep1[i][j]=vep1[i][k];
                        vep1[i][k]=p;
                    }
                }
            }
        }
        for(int i=0;i<n;i++){
            for(int j=0;j<n;j++){
                System.out.printf("\t%f\t",vep1[i][j]);
            }
            System.out.println("\n");
        }
    }
    static void individus(double[][] matR,double[][] vep1,double[][] matCOmP){
        System.out.println("\n\n\t ************************* individus************************* \n\n");
        for(int i=0;i<m;i++){
            for(int j=0;j<n;j++){
                matCOmP[i][j]=0;
                for(int k=0;k<n;k++){
                    matCOmP[i][j]=matCOmP[i][j]+(matR[i][k]*vep1[j][k]);
                }
            }
        }
        for(int i=0;i<m;i++){
            for(int j=0;j<n;j++){
                System.out.printf("\t%f\t",matCOmP[i][j]);
            }
            System.out.println("\n");
        }       
    }
    static void valeursPr(double[][] v1,double[][] corr){
        System.out.println("\n\n\t************************ Les valeurs Propre ************************* \n\n");
        Matrix m1= new Matrix(corr);
        double[][] v2=new double[n][n];
        v2= m1.eig().getD().getArray();
        for(int i = 0; i <n; i++){
            v1[i][i]=v2[i][i];
            System.out.printf("\t%f\t",v1[i][i]);
        }
    }
    static void Coefficients_correlation(double[][] v1,double[][] matC,double[][] vep2){
        System.out.println("\n\n\t************************* Coefficients de corrélation ************************* \n\n");
        for(int i=0; i< n ; i++){
            for(int j=0;j<n;j++){
                matC[i][j]=Math.sqrt(v1[i][i])*vep2[j][i];
            }
        }
        for(int i=0;i<n;i++){
            for(int j=0;j<n;j++){
                System.out.printf("\t%f\t",matC[i][j]);
            }
            System.out.println("\n");
        }
    }
    


    public static void main(String[] args) throws Exception {
        double[][] matNT = new double[n][m];
        double[][] matC = new double[m][n];
        double[][] matR = new double[m][n];
        double[][] corr = new double[n][n];
        double[][] matN = new double[m][n];
        double[][] matINd=new double[m][n];
        double [][]vep1=new double[n][n];
        double [][]vep2=new double[n][n];
        double[][] v1 = new double[n][n];
        double[] V1=new double[15];
        double[] V2 = new double[15];
        double[] V3= new double[15];
        //double[] x= new double[n];
        //double[] y= new double[n];
        //double[] z= new double[n];
        double[][] tab = {
            {5.1f, 3.5f, 1.4f},
            {4.9f, 3.0f, 1.4f},
            {4.7f, 3.2f, 1.3f},
            {4.6f, 3.1f, 1.5f},
            {5.0f, 3.6f, 1.4f},
            {7.0f, 3.2f, 4.7f},
            {6.4f, 3.2f, 4.5f},
            {6.9f, 3.1f, 4.9f},
            {5.5f, 2.3f, 4.0f},
            {6.5f, 2.8f, 4.6f},
            {6.3f, 3.3f, 6.0f},
            {5.8f, 2.7f, 5.1f},
            {7.1f, 3.0f, 5.9f},
            {6.3f, 2.9f, 5.6f},
            {6.5f, 3.0f,5.8f}
        };
        affichermat(tab);
        centrée(tab,matC);
        reduit(matC, matR);
        Normees(matR, matN);
        transpose(matN, matNT);
        CORRELATION(matNT, matN, corr);
        vecteurPr(corr,vep2);
        vecteurPr2(vep1, vep2);
        vecteurPrTr(vep1);
        individus(matR, vep1, matINd);
        valeursPr(v1, corr);
        Coefficients_correlation(v1, matC, vep2);
        MatToVec(matINd,V1,V2,V3);
        g3D(V1, V2, V3);
        
    }
}
