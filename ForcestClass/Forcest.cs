using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ForcestClass
{
    public class Forcest
    {
        Random rand_seed =null;
        int rnd = 0;
        public Forcest()
        {

            int max = 1000;
            int rnd = int.MinValue;
            decimal _base = (decimal)long.MaxValue;
            byte[] rndSeries = new byte[8];
            System.Security.Cryptography.RNGCryptoServiceProvider rng
                = new System.Security.Cryptography.RNGCryptoServiceProvider();
            rng.GetBytes(rndSeries);
            //不含100需去掉+1 
            rnd = (int)(Math.Abs(BitConverter.ToInt64(rndSeries, 0)) / _base * (max + 1));

            rand_seed = new Random(rnd);
        }
        private void MakeRandMat(ref Queue<Queue<double>> mat,int row ,int col)
        {
            Random rand = new Random(rand_seed.Next());
            for (int i = 0; i < mat.Count; i++)
            {
                if (i > row - 1)
                    mat.Dequeue();
                else
                {
                    Queue<double> tempQueue = new Queue<double>();
                    for (int j = 0; j < col; j++)
                    {
                        tempQueue.Enqueue(rand.Next(0, 1000)/1000.0);
                    }
                    mat.Enqueue(tempQueue);
                }
            }
            //for (int i = 0; i < ; i++)
            //{
            //    for (int j = 0; j < mat.ElementAt(i).Count; j++)
            //    {
            //        mat.ElementAt(i). = rand.Next(0,1000);
            //    }
                
            //}
            
            //double u1, u2, v1 = 0, v2 = 0, s = 0, z1 = 0, z2 = 0;
            //while (s > 1 || s == 0)
            //{
            //    u1 = rand.NextDouble();
            //    u2 = rand.NextDouble();
            //    v1 = 2 * u1 - 1;
            //    v2 = 2 * u2 - 1;
            //    s = v1 * v1 + v2 * v2;
            //}
            //z1 = Math.Sqrt(-2 * Math.Log(s) / s) * v1;
            //z2 = Math.Sqrt(-2 * Math.Log(s) / s) * v2;
            //y = new double[] { z1, z2 };
            //return y; //返回两个服从正态分布N(0,1)的随机数z0 和 z1
        }
        ///////////////////////////////////////////////////////////////////////


        //////////预测算法///////////////////////////////////////////////////////
        //vector< vector< double > > Input;////////////训练输入样本
        //vector< vector< double > > Output;//////////训练输出样本
        //vector< vector< double > > PreInput;////////预测输入样本
        //vector< vector< double > > PreOutput; ///////预测输出结果
        //DataTable Input = new DataTable();
        //DataTable Output = new DataTable();
        //DataTable PreInput = new DataTable();
        //DataTable PreOutput = new DataTable();
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        /////本人为动态改变预测体，采用vector向量结构计算，可随时改变数组结构///////////////////////////
        /////////vector的使用需添加标准库头文件<vertor>，并定义命名空间，using namespace std;///
        ///////////////规避double数组只能静态定义的特性///////////////////////////////////////////////////////////////////////////////////
        ///////////调试中发现double二维及以上常规数组无法进行指针与数组名之间的转化/////////////////////
        /////////////具体错误cannot convert double** to double[][]//////////////////////////////////////// 
        void Forecast()
        {


            List<List<double>> Input = new List<List<double>>();
            List<List<double>> Output = new List<List<double>>();
            List<List<double>> PreInput = new List<List<double>>();
            List<List<double>> PreOutput = new List<List<double>>();
            int input_col;//定义训练输入样本列，数据计算以列为单位，与MATLAB中列向量同意义
            int input_row;//定义训练输入样本行
            int inputtest_row;//定义训练输出样本行
            int output_col;//定义预测输入样本列，即需要预测的数据，1列
            int output_row;//定义预测输入样本行
            int outputtest_row;//定义预测输出样本行，1列
            int maxgen; //迭代次数
            int n; //隐形节点个数	
            //double[,] Input = new double[20,1];
            ///////////////以下为个训练向量与预测向量结构初始化，预测时需进行数据有效化填充//////////////////////////////// 
            int i = 0;
            //Input = new List<double>(input_row);
            //Input = new List<List<double>>;
            Input = new List<List<double>>(input_row);
            for (i = 0; i < input_row; i++)
            {
                Input[i] = new List<double>(input_col);
            }
            Output = new List<List<double>>(inputtest_row);
            for (i = 0; i < inputtest_row; i++)
            {
                Output[i] = new List<double>(input_col);
            }
            PreInput = new List<List<double>>(output_row);
            for (i = 0; i < output_row; i++)
            {
                PreInput[i] = new List<double>(output_col);
            }
            PreOutput = new List<List<double>>(outputtest_row);
            for (i = 0; i < outputtest_row; i++)
            {
                PreOutput[i] = new List<double>(output_col);
            }



            int M = input_col;//输入节点个数
            int N = output_col;//输出节点个数
            double lr1 = 0.01; ////学习概率
            double lr2 = 0.001; ////学习概率
            int j, k, kk, kkk;
            //srand((unsigned)time( NULL ));
            double value = 0.0;
            string str = "";
            List<List<double>> Wjk = new List<List<double>>(n);
            for (i = 0; i < n; i++)
            {
                Wjk[i] = new List<double>(input_col);
            }
            List<List<double>> Wjk_1 = new List<List<double>>(n);
            for (i = 0; i < n; i++)
            {
                Wjk_1[i] = new List<double>(input_col);
            }
            List<List<double>> Wjk_2 = new List<List<double>>(n);
            for (i = 0; i < n; i++)
            {
                Wjk_2[i] = new List<double>(input_col);
            }

            for (i = 0; i < n; i++)
            {
                for (j = 0; j < input_col; j++)
                {
                    Wjk[i][j] = randn(2);
                    Wjk_1[i][j] = Wjk[i][j];
                    Wjk_2[i][j] = Wjk_1[i][j];
                }
            }
            List<List<double>> Wij = new List<List<double>>(output_col);
            for (i = 0; i < output_col; i++)
            {
                Wij[i] = new List<double>(n);
            }
            List<List<double>> Wij_1 = new List<List<double>>(output_col);
            for (i = 0; i < output_col; i++)
            {
                Wij_1[i] = new List<double>(n);
            }
            List<List<double>> Wij_2 = new List<List<double>>(output_col);
            for (i = 0; i < output_col; i++)
            {
                Wij_2[i] = new List<double>(n);
            }


            for (i = 0; i < output_col; i++)
            {
                for (j = 0; j < n; j++)
                {
                    Wij[i][j] = randn(2);
                    Wij_1[i][j] = Wij[i][j];
                    Wij_2[i][j] = Wij_1[i][j];
                }
            }

            List<List<double>> a = new List<List<double>>(1);
            for (i = 0; i < 1; i++)
            {
                a[i] = new List<double>(n);
            }
            List<List<double>> a_1 = new List<List<double>>(1);
            for (i = 0; i < 1; i++)
            {
                a_1[i] = new List<double>(n);
            }
            List<List<double>> a_2 = new List<List<double>>(1);
            for (i = 0; i < 1; i++)
            {
                a_2[i] = new List<double>(n);
            }

            for (i = 0; i < 1; i++)
            {
                for (j = 0; j < n; j++)
                {
                    a[i][j] = randn(2);
                    a_1[i][j] = a[i][j];
                    a_2[i][j] = a_1[i][j];
                }
            }

            List<List<double>> b = new List<List<double>>(1);
            for (i = 0; i < 1; i++)
            {
                b[i] = new List<double>(n);
            }
            List<List<double>> b_1 = new List<List<double>>(1);
            for (i = 0; i < 1; i++)
            {
                b_1[i] = new List<double>(n);
            }
            List<List<double>> b_2 = new List<List<double>>(1);
            for (i = 0; i < 1; i++)
            {
                b_2[i] = new List<double>(n);
            }
            for (i = 0; i < 1; i++)
            {
                for (j = 0; j < n; j++)
                {
                    b[i][j] = randn(2);
                    b_1[i][j] = b[i][j];
                    b_2[i][j] = b_1[i][j];
                }
            }
            List<List<double>> y = new List<List<double>>(1);
            for (i = 0; i < 1; i++)
            {
                y[i] = new List<double>(output_col);
            }
            List<List<double>> net = new List<List<double>>(1);
            for (i = 0; i < 1; i++)
            {
                net[i] = new List<double>(n);
            }
            List<List<double>> net_ab = new List<List<double>>(1);
            for (i = 0; i < 1; i++)
            {
                net_ab[i] = new List<double>(n);
            }
            List<List<double>> d_Wjk = new List<List<double>>(n);
            for (i = 0; i < n; i++)
            {
                d_Wjk[i] = new List<double>(input_col);
            }
            List<List<double>> d_Wij = new List<List<double>>(output_col);
            for (i = 0; i < output_col; i++)
            {
                d_Wij[i] = new List<double>(n);
            }
            List<List<double>> d_a = new List<List<double>>(1);
            for (i = 0; i < 1; i++)
            {
                d_a[i] = new List<double>(n);
            }
            List<List<double>> d_b = new List<List<double>>(1);
            for (i = 0; i < 1; i++)
            {
                d_b[i] = new List<double>(n);
            }
            ////// 输入输出数据归一化
            double tempmax = 0.0, tempmin = 0.0;
            List<double> input_max = new List<double>(input_col);
            List<double> input_min = new List<double>(input_col);
            List<List<double>> Input1 = new List<List<double>>(input_row);
            for (i = 0; i < input_row; i++)
            {
                Input1[i] = new List<double>(input_col);
            }
            for (i = 0; i < input_col; i++)
            {
                tempmax = Input[0][i];
                tempmin = Input[0][i];
                for (j = 0; j < input_row; j++)
                {
                    if (Input[j][i] > tempmax)
                        tempmax = Input[j][i];
                    if (Input[j][i] < tempmin)
                        tempmin = Input[j][i];
                }
                input_max[i] = tempmax;
                input_min[i] = tempmin;
                for (j = 0; j < input_row; j++)
                {
                    if (input_max[i] == input_min[i])
                        Input1[j][i] = 0;
                    else
                        Input1[j][i] = 2 * (Input[j][i] - input_min[i]) / (input_max[i] - input_min[i]) - 1;
                }
            }
            List<double> output_max = new List<double>(output_col);
            List<double> output_min = new List<double>(output_col);
            List<List<double>> PreInput1 = new List<List<double>>(output_row);
            for (i = 0; i < output_row; i++)
            {
                PreInput1[i] = new List<double>(output_col);
            }
            for (i = 0; i < output_col; i++)
            {
                tempmax = PreInput[0][i];
                tempmin = PreInput[0][i];
                for (j = 0; j < output_row; j++)
                {
                    if (PreInput[j][i] > tempmax)
                        tempmax = PreInput[j][i];
                    if (PreInput[j][i] < tempmin)
                        tempmin = PreInput[j][i];
                }
                output_max[i] = tempmax;
                output_min[i] = tempmin;
                for (j = 0; j < output_row; j++)
                {
                    if (output_max[i] == output_min[i])
                        PreInput1[j][i] = 0;
                    else
                        PreInput1[j][i] = 2 * (PreInput[j][i] - output_min[i]) / (output_max[i] - output_min[i]) - 1;
                }
            }

            ////// 网络训练
            List<double> x = new List<double>(input_col);
            List<double> yqw = new List<double>(output_col);
            List<double> error = new List<double>(maxgen);
            double temp = 0.0;
            for (i = 0; i < maxgen; i++)
            {
                ////误差累计
                error[i] = 0.0;
                for (kk = 0; kk < input_row; kk++)
                {
                    for (kkk = 0; kkk < input_col; kkk++)
                    {
                        x[kkk] = Input1[kk][kkk];
                    }
                    for (kkk = 0; kkk < output_col; kkk++)
                    {
                        yqw[kkk] = PreInput1[kk][kkk];
                    }
                    for (j = 0; j < n; j++)
                    {
                        for (k = 0; k < input_col; k++)
                        {
                            net[0][j] = net[0][j] + Wjk[j][k] * x[k];
                            net_ab[0][j] = (net[0][j] - b[0][j]) / a[0][j];
                        }
                        temp = mymorlet(net_ab[0][j]);
                        for (k = 0; k < output_col; k++)
                        {
                            y[0][k] = y[0][k] + Wij[k][j] * temp;
                        }
                    }
                    for (j = 0; j < output_col; j++)
                    {
                        temp = temp + abs(yqw[j] - y[0][j]);
                    }
                    error[i] = error[i] + temp;
                    for (j = 0; j < n; j++)
                    {
                        temp = mymorlet(net_ab[0][j]);
                        for (k = 0; k < output_col; k++)
                        {
                            d_Wij[k][j] = d_Wij[k][j] - (yqw[k] - y[0][k]) * temp;
                        }
                        temp = d_mymorlet(net_ab[0][j]);
                        for (k = 0; k < input_col; k++)
                        {
                            for (kkk = 0; kkk < output_col; kkk++)
                            {
                                d_Wjk[j][k] = d_Wjk[j][k] + (yqw[kkk] - y[0][kkk]) * Wij[kkk][j];
                            }
                            d_Wjk[j][k] = -d_Wjk[j][k] * temp * x[k] / a[0][j];
                        }
                        for (k = 0; k < output_col; k++)
                        {
                            d_b[0][j] = d_b[0][j] + (yqw[k] - y[0][k]) * Wij[k][j];
                        }
                        d_b[0][j] = d_b[0][j] * temp / a[0][j];
                        for (k = 0; k < output_col; k++)
                        {
                            d_a[0][j] = d_a[0][j] + (yqw[k] - y[0][k]) * Wij[k][j];
                        }
                        d_a[0][j] = d_a[0][j] * temp * ((net[0][j] - b[0][j]) / b[0][j]) / a[0][j];
                    }

                    ///权值参数更新
                    for (j = 0; j < n; j++)
                    {
                        for (k = 0; k < input_col; k++)
                        {
                            Wjk[j][k] = Wjk[j][k] - lr1 * d_Wjk[j][k];
                            d_Wjk[j][k] = 0.0;
                            Wjk_1[j][k] = Wjk[j][k];
                            Wjk_2[j][k] = Wjk_1[j][k];
                        }
                        b[0][j] = b[0][j] - lr2 * d_b[0][j];
                        a[0][j] = a[0][j] - lr2 * d_a[0][j];
                        d_a[0][j] = 0.0;
                        d_b[0][j] = 0.0;
                        net[0][j] = 0.0;
                        net_ab[0][j] = 0.0;
                        a_1[0][j] = a[0][j];
                        a_2[0][j] = a_1[0][j];
                        b_1[0][j] = b[0][j];
                        b_1[0][j] = b_1[0][j];
                    }
                    for (j = 0; j < output_col; j++)
                    {
                        for (k = 0; k < n; k++)
                        {
                            Wij[j][k] = Wij[j][k] - lr1 * d_Wij[j][k];
                            d_Wij[j][k] = 0.0;
                            Wij_1[j][k] = Wij[j][k];
                            Wij_2[j][k] = Wij_1[j][k];
                        }
                        y[0][j] = 0.0;
                    }
                }
            }


            ////// 网络预测
            //////预测输入归一化
            List<double> inputtest_max = new List<double>(input_col);
            List<double> inputtest_min = new List<double>(input_col);
            List<List<double>> OutputAsInput1 = new List<List<double>>(inputtest_row);
            for (i = 0; i < inputtest_row; i++)
            {
                OutputAsInput1[i] = new List<double>(input_col);
            }
            //double inputtest_min[input_col] = {{0.0}};
            for (i = 0; i < input_col; i++)
            {
                tempmax = Output[0][i];
                tempmin = Output[0][i];
                for (j = 0; j < inputtest_row; j++)
                {
                    if (Output[j][i] > tempmax)
                        tempmax = Output[j][i];
                    if (Output[j][i] < tempmin)
                        tempmin = Output[j][i];
                }
                inputtest_max[i] = tempmax;
                inputtest_min[i] = tempmin;
                for (j = 0; j < inputtest_row; j++)
                {
                    if (inputtest_max[i] == inputtest_min[i])
                        OutputAsInput1[j][i] = 0.5 * (input_max[i] + input_min[i]);
                    else
                        OutputAsInput1[j][i] = (Output[j][i] - inputtest_min[i]) / (inputtest_max[i] - inputtest_min[i]) * (input_max[i] - input_min[i]) + input_min[i];
                }
            }

            /////网络预测
            List<double> x_test = new List<double>(input_col);
            List<List<double>> yuce = new List<List<double>>(1);
            for (i = 0; i < 1; i++)
            {
                yuce[i] = new List<double>(inputtest_row);
            }
            for (i = 0; i < inputtest_row; i++)
            {
                for (j = 0; j < input_col; j++)
                {
                    x_test[j] = OutputAsInput1[i][j];
                }
                for (j = 0; j < n; j++)
                {
                    for (k = 0; k < input_col; k++)
                    {
                        net[0][j] = net[0][j] + Wjk[j][k] * x_test[k];
                        net_ab[0][j] = (net[0][j] - b[0][j]) / a[0][j];
                    }
                    temp = mymorlet(net_ab[0][j]);
                    for (k = 0; k < output_col; k++)
                    {
                        y[0][k] = y[0][k] + Wij[k][j] * temp;
                    }
                }
                yuce[0][i] = y[0][k - 1];
                for (j = 0; j < output_col; j++)
                {
                    y[0][j] = 0.0;
                }
                for (j = 0; j < n; j++)
                {
                    net[0][j] = 0.0;
                    net_ab[0][j] = 0.0;
                }
            }

            List<double> yun = new List<double>(outputtest_row);
            for (i = 0; i < outputtest_row; i++)
            {
                yun[i] = (yuce[0][i] + 1) / 2 * (output_max[0] - output_min[0]) + output_min[0];
                PreOutput[i][output_col - 1] = yun[i];
            }
            return;
        }


        private void forcest()
        {
        }
    }
}
