using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Data;

namespace ClassForecast
{
    public class Forecast
    {
        Log logtxt = null;
        private static Random rand_seed = null;
        private Queue<Queue<double>> input = null;//各个序列以列为Queue单位
        private Queue<Queue<double>> input_test = null;
        private Queue<double> output = null;
        //private Queue<double> output_test = null;
        private Queue<double> source_data = null;
        private int items_count = new int();
        public bool IsFilled = false;

        public int input_row = new int();
        public int input_test_row = new int();
        public int input_col = new int();
        public int sample_num = new int();//sample_row应该等于input_row
        //private int forecast_num = new int();


        public Forecast()
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
            logtxt = new Log(AppDomain.CurrentDomain.BaseDirectory + @"/log/MatlabLog.txt");
            input = new Queue<Queue<double>>();
            input_test = new Queue<Queue<double>>();
            output = new Queue<double>();
            //output_test = new Queue<double>();
            source_data = new Queue<double>();
            items_count = 0;
            input_row = 0;
            input_test_row = 0;
            input_col = 0;
            sample_num = 0;
            //forecast_num = 0;
        }

        private static double[,] MakeRandDbMat(int row, int col, double min = 0.0, double max = 1.0)
        {
            Random rand = new Random(rand_seed.Next());
            double[,] result = new double[row, col];
            for (int i = 0; i < row; i++)
            {
                for (int j = 0; j < col; j++)
                {
                    result[i, j] = rand.Next((int)min * 1000, (int)max * 1000) / 1000.0;
                }
            }
            return result;
        }

        private static double[,] GetDbMatMaxMin(Queue<Queue<double>> mat)//获取矩阵列最大值最小值，最大在前，最小在后
        {
            double[,] result = new double[2, mat.Count];
            for (int i = 0; i < mat.Count; i++)
            {
                result[0, i] = mat.ElementAt(i).Max();
                result[1, i] = mat.ElementAt(i).Min();
            }
            return result;
        }

        private static double[] GetDbMatMaxMin(Queue<double> mat)//获取矩阵列最大值最小值，最大在前，最小在后
        {
            double[] result = new double[2];
            result[0] = mat.Max();
            result[1] = mat.Min();
            return result;
        }

        private static double[,] DbMapMinMax(Queue<double> mat, double[] maxmin)//二维矩阵归一化,返回数组
        {
            double[,] result = new double[mat.Count, 1];
            for (int i = 0; i < mat.Count; i++)
            {
                for (int j = 0; j < 1; j++)
                {
                    if (maxmin[0] == maxmin[1])
                        result[i, j] = 0.0;
                    else
                        result[i, j] = ((mat.ElementAt(i) - maxmin[1]) * 2) / (maxmin[0] - maxmin[1]) - 1;
                }
            }
            return result;
        }

        private static double[,] DbMapMinMax(Queue<Queue<double>> mat, double[,] maxmin)//二维矩阵归一化,返回数组
        {
            double[,] result = new double[mat.ElementAt(0).Count, mat.Count];
            for (int i = 0; i < mat.Count; i++)
            {
                for (int j = 0; j < mat.ElementAt(0).Count; j++)
                {
                    if (maxmin[0, i] == maxmin[1, i])
                        result[j, i] = 0.0;
                    else
                        result[j, i] = ((mat.ElementAt(i).ElementAt(j) - maxmin[1, i]) * 2) / (maxmin[0, i] - maxmin[1, i]) - 1;
                }
            }
            return result;
        }

        private static double[,] ReserveMaxMin(double[,] mat, double[,] maxmin)
        {
            int row = mat.Length / maxmin.Length * 2;
            int col = maxmin.Length / 2;
            double[,] result = new double[row, col];
            for (int i = 0; i < row; i++)
            {
                for (int j = 0; j < col; j++)
                {
                    result[i, j] = (mat[i, j] + 1) * (maxmin[0, j] - maxmin[1, j]) / 2 + maxmin[1, j];
                }
            }
            return result;
        }

        private static double[] ReserveMaxMin(double[] mat, double[] maxmin)
        {
            int row = mat.Length;
            double[] result = new double[row];
            for (int i = 0; i < row; i++)
            {
                result[i] = (mat[i] + 1) * (maxmin[0] - maxmin[1]) / 2 + maxmin[1];
            }
            return result;
        }

        private static double[,] ApplyMapMinMax(Queue<Queue<double>> mat, double[,] maxmin)
        {
            double[,] result = new double[mat.ElementAt(0).Count, mat.Count];
            for (int i = 0; i < mat.Count; i++)
            {
                for (int j = 0; j < mat.ElementAt(0).Count; j++)
                {
                    if (maxmin[0, i] == maxmin[1, i])
                        result[j, i] = 0.0;
                    else
                        result[j, i] = ((mat.ElementAt(i).ElementAt(j) - maxmin[1, i]) * 2) / (maxmin[0, i] - maxmin[1, i]) - 1;
                }
            }
            return result;
        }

        private static double mymorlet(double data)
        {
            return Math.Exp(-Math.Pow(data, 2) / 2) * Math.Cos(1.75 * data);
        }

        private static double d_mylorlet(double data)
        {
            return -1.75 * Math.Sin(1.75 * data) * Math.Exp(-Math.Pow(data, 2) / 2) - data * Math.Cos(1.75 * data) * Math.Exp(-Math.Pow(data, 2) / 2);
        }

        public bool AddDataToSource(double data)
        {
            AddDataToQueue(data);
            source_data.Enqueue(data);
            items_count++;
            if (items_count > 500)
            {
                source_data.Dequeue();
                items_count--;
            }
            if (items_count >= (input_row + input_test_row) * input_col)
            {
                IsFilled = true;
                return true;
            }
            else
            {
                IsFilled = false;
                return false;
            }
        }

        private void AddDataToQueue(double data)
        {
            while (input.Count < input_col)
            {
                Queue<double> tempQueue = new Queue<double>();
                tempQueue.Enqueue(0.0);
                input.Enqueue(tempQueue);
            }
            while (input_test.Count < input_col)
            {
                Queue<double> tempQueue = new Queue<double>();
                tempQueue.Enqueue(0.0);
                input_test.Enqueue(tempQueue);
            }
            if (input_row != input.ElementAt(0).Count || input_test_row != input_test.ElementAt(0).Count)
            {
                if (items_count >= (input_row + input_test_row) * input_col)
                {
                    input.Clear();
                    input_test.Clear();
                    output.Clear();
                    int i = 0;
                    int FirstPos = items_count - (input_row + input_test_row) * input_col;
                    Queue<double> tempQueue = new Queue<double>();
                    for (int j = 0; j < input_col; j++)
                    {
                        tempQueue.Clear();
                        for (int k = 0; k < input_row; k++)
                        {
                            tempQueue.Enqueue(source_data.ElementAt(FirstPos + i++));
                        }
                        input.Enqueue(tempQueue);
                        tempQueue.Clear();
                        for (int k = 0; k < input_test_row; k++)
                        {
                            tempQueue.Enqueue(source_data.ElementAt(FirstPos + i++));
                        }
                        input_test.Enqueue(tempQueue);
                    }
                }
                else
                {
                    input.Clear();
                    input_test.Clear();
                    output.Clear();
                    int i = 0, l = 0;
                    int ZeroCount = (input_row + input_test_row) * input_col - items_count;
                    Queue<double> tempQueue = new Queue<double>();
                    for (int j = 0; j < input_col; j++)
                    {
                        tempQueue.Clear();
                        for (int k = 0; k < input_row; k++)
                        {
                            if (l < ZeroCount)
                            {
                                tempQueue.Enqueue(new double());
                                l++;
                            }
                            else
                                tempQueue.Enqueue(source_data.ElementAt(i++));
                        }
                        input.Enqueue(tempQueue);
                        tempQueue.Clear();
                        for (int k = 0; k < input_test_row; k++)
                        {
                            if (l < ZeroCount)
                            {
                                tempQueue.Enqueue(new double());
                                l++;
                            }
                            else
                                tempQueue.Enqueue(source_data.ElementAt(i++));
                        }
                        input_test.Enqueue(tempQueue);
                    }

                }
                output.Clear();
                for (int i = 0; i < sample_num; i++)
                {
                    if (items_count >= sample_num)
                        output.Enqueue(source_data.ElementAt(items_count - sample_num + i));
                    else
                    {
                        if (i < sample_num - items_count)
                        {
                            output.Enqueue(new double());
                        }
                        else
                            output.Enqueue(source_data.ElementAt(items_count - sample_num + i));
                    }
                }
            }

            for (int j = 0; j < input_col - 1; j++)
            {
                input.ElementAt(j).Dequeue();
                input.ElementAt(j).Enqueue(input_test.ElementAt(j).Peek());
                input_test.ElementAt(j).Dequeue();
                input_test.ElementAt(j).Enqueue(input.ElementAt(j + 1).Peek());
            }
            input.ElementAt(input_col - 1).Dequeue();
            input.ElementAt(input_col - 1).Enqueue(input_test.ElementAt(input_col - 1).Peek());
            input_test.ElementAt(input_col - 1).Dequeue();
            input_test.ElementAt(input_col - 1).Enqueue(data);

            output.Enqueue(data);
            while (output.Count > sample_num)
                output.Dequeue();
        }

        public double[] MakeForecast()
        {
            //            
            //% %% 清空环境变量
            //% clc
            //% clear
            //% 
            //% %% 网络参数配置
            //% load traffic_flux input output input_test output_test
            //function wavenn(input,output,input_test,output_test)

            //M=size(input,2); %输入节点个数
            //N=size(output,2); %输出节点个数
            if (!IsFilled)
            {
                return null;
            }
            int M = input.Count;
            int N = 1;//output.Count

            int n = 6;
            double lr1 = 0.01;
            double lr2 = 0.001;
            int maxgen = 100;

            //n=6; %隐形节点个数
            //lr1=0.01; %学习概率
            //lr2=0.001; %学习概率
            //maxgen=100; %迭代次数


            double[,] Wjk = MakeRandDbMat(n, M);
            double[,] Wjk_1 = Wjk;
            double[,] Wjk_2 = Wjk;

            double[,] Wij = MakeRandDbMat(N, n);
            double[,] Wij_1 = Wjk;
            double[,] Wij_2 = Wjk;

            double[,] a = MakeRandDbMat(1, n);
            double[,] a_1 = a;
            double[,] a_2 = a;

            double[,] b = MakeRandDbMat(1, n);
            double[,] b_1 = b;
            double[,] b_2 = b;
            //%权值初始化
            //Wjk=randn(n,M);Wjk_1=Wjk;Wjk_2=Wjk_1;
            //Wij=randn(N,n);Wij_1=Wij;Wij_2=Wij_1;
            //a=randn(1,n);a_1=a;a_2=a_1;
            //b=randn(1,n);b_1=b;b_2=b_1;
            double y = 0.0;
            double[,] net = MakeRandDbMat(1, n, 0, 0);
            double[,] net_ab = MakeRandDbMat(1, n, 0, 0);
            //%节点初始化
            //y=zeros(1,N);
            //net=zeros(1,n);
            //net_ab=zeros(1,n);
            double[,] d_Wjk = MakeRandDbMat(n, M, 0, 0);
            double[,] d_Wij = MakeRandDbMat(N, n, 0, 0);
            double[,] d_a = MakeRandDbMat(1, n, 0, 0);
            double[,] d_b = MakeRandDbMat(1, n, 0, 0);
            //%权值学习增量初始化
            //d_Wjk=zeros(n,M);
            //d_Wij=zeros(N,n);
            //d_a=zeros(1,n);
            //d_b=zeros(1,n);
            double[,] inputps = GetDbMatMaxMin(input);
            double[,] inputn = DbMapMinMax(input, inputps);
            double[] outputps = GetDbMatMaxMin(output);
            double[,] outputn = DbMapMinMax(output, outputps);
            //%% 输入输出数据归一化


            //[inputn,inputps]=mapminmax(input');
            //[outputn,outputps]=mapminmax(output'); 
            //inputn=inputn';
            //outputn=outputn';
            double[] error = new double[maxgen];
            for (int i = 0; i < maxgen; i++)
            {
                error[i] = 0;
                for (int kk = 0; kk < input.ElementAt(0).Count; kk++)
                {
                    double[] x = new double[input.Count];
                    double yqw = new double();
                    for (int j = 0; j < input.Count; j++)
                    {
                        x[j] = inputn[kk, j];
                    }
                    yqw = outputn[kk, 0];
                    for (int j = 0; j < n; j++)
                    {
                        for (int k = 0; k < M; k++)
                        {
                            net[0, j] = net[0, j] + Wjk[j, k] * x[k];
                            net_ab[0, j] = (net_ab[0, j] - b[0, j]) / a[0, j];
                        }
                        double temp = mymorlet(net_ab[0, j]);
                        for (int k = 0; k < N; k++)
                        {
                            y = y + Wij[k, j] * temp;
                        }
                    }
                    error[i] = error[i] + Math.Abs(yqw - y);

                    for (int j = 0; j < n; j++)
                    {
                        double temp = mymorlet(net_ab[0, j]);
                        for (int k = 0; k < N; k++)
                        {
                            d_Wij[k, j] = d_Wij[k, j] - (yqw - y) * temp;
                        }
                        temp = d_mylorlet(net_ab[0, j]);
                        for (int k = 0; k < M; k++)
                        {
                            for (int l = 0; l < N; l++)
                            {
                                d_Wjk[j, k] = d_Wjk[j, k] + (yqw - y) * Wij[l, j];
                            }
                            d_Wjk[j, k] = -d_Wjk[j, k] * temp * x[k] / a[0, j];
                        }
                        d_b[0, j] = d_b[0, j] + (yqw - y) * Wij[0, j];
                        d_b[0, j] = d_b[0, j] * temp / a[0, j];
                        d_a[0, j] = d_a[0, j] + (yqw - y) * Wij[0, j];
                        d_a[0, j] = d_a[0, j] * temp * ((net[0, j] - b[0, j]) / b[0, j]) / a[0, j];
                    }

                    //权值参数更新
                    for (int j = 0; j < N; j++)
                    {
                        for (int k = 0; k < n; k++)
                        {
                            Wij[j, k] = Wij[j, k] - lr1 * d_Wij[j, k];
                        }
                    }
                    for (int j = 0; j < n; j++)
                    {
                        for (int k = 0; k < M; k++)
                        {
                            Wjk[j, k] = Wjk[j, k] - lr1 * d_Wjk[j, k];
                        }
                    }
                    for (int j = 0; j < 1; j++)
                    {
                        for (int k = 0; k < n; k++)
                        {
                            b[j, k] = b[j, k] - lr2 * d_b[j, k];
                        }
                    }
                    for (int j = 0; j < 1; j++)
                    {
                        for (int k = 0; k < n; k++)
                        {
                            a[j, k] = a[j, k] - lr2 * d_a[j, k];
                        }
                    }
                    d_Wjk = MakeRandDbMat(n, M, 0, 0);
                    d_Wij = MakeRandDbMat(N, n, 0, 0);
                    d_a = MakeRandDbMat(1, n, 0, 0);
                    d_b = MakeRandDbMat(1, n, 0, 0);

                    y = 0;
                    net = MakeRandDbMat(1, n, 0, 0);
                    net_ab = MakeRandDbMat(1, n, 0, 0);

                    Wjk_1 = Wjk;
                    Wjk_2 = Wjk;
                    Wij_1 = Wij;
                    Wij_2 = Wij;
                    a_1 = a;
                    a_2 = a;
                    b_1 = b;
                    b_2 = b;
                }
            }
            double[,] xx = DbMapMinMax(input_test, inputps);
            double[] yuce = new double[input_test.ElementAt(0).Count];
            for (int i = 0; i < input_test.ElementAt(0).Count; i++)
            {
                double[] xx_test = new double[input_test.ElementAt(0).Count];
                for (int j = 0; j < input_test.Count; j++)
                {
                    xx_test[j] = input_test.ElementAt(j).ElementAt(i);
                }
                for (int j = 0; j < n; j++)
                {
                    for (int k = 0; k < M; k++)
                    {
                        net_ab[0, j] = net[0, j] + Wjk[j, k] * xx_test[k];
                        net_ab[0, j] = (net_ab[0, j] - b[0, j]) / a[0, j];
                    }
                    double temp = mymorlet(net_ab[0, j]);
                    for (int k = 0; k < N; k++)
                    {
                        y = y + Wij[k, j] * temp;
                    }
                }
                yuce[i] = y;
                y = 0.0;
                net = MakeRandDbMat(1, n, 0, 0);
                net_ab = MakeRandDbMat(1, n, 0, 0);
            }
            double[] yun = ReserveMaxMin(yuce, outputps);
            return yun;
        }

        public DataTable MakeForecast(int glag)
        {
            //            
            //% %% 清空环境变量
            //% clc
            //% clear
            //% 
            //% %% 网络参数配置
            //% load traffic_flux input output input_test output_test
            //function wavenn(input,output,input_test,output_test)

            //M=size(input,2); %输入节点个数
            //N=size(output,2); %输出节点个数
            if (!IsFilled)
            {
                return null;
            }
            int M = input.Count;
            int N = 1;//output.Count

            int n = 6;
            double lr1 = 0.01;
            double lr2 = 0.001;
            int maxgen = 100;

            //n=6; %隐形节点个数
            //lr1=0.01; %学习概率
            //lr2=0.001; %学习概率
            //maxgen=100; %迭代次数


            double[,] Wjk = MakeRandDbMat(n, M);
            double[,] Wjk_1 = Wjk;
            double[,] Wjk_2 = Wjk;

            double[,] Wij = MakeRandDbMat(N, n);
            double[,] Wij_1 = Wjk;
            double[,] Wij_2 = Wjk;

            double[,] a = MakeRandDbMat(1, n);
            double[,] a_1 = a;
            double[,] a_2 = a;

            double[,] b = MakeRandDbMat(1, n);
            double[,] b_1 = b;
            double[,] b_2 = b;
            //%权值初始化
            //Wjk=randn(n,M);Wjk_1=Wjk;Wjk_2=Wjk_1;
            //Wij=randn(N,n);Wij_1=Wij;Wij_2=Wij_1;
            //a=randn(1,n);a_1=a;a_2=a_1;
            //b=randn(1,n);b_1=b;b_2=b_1;
            double y = 0.0;
            double[,] net = MakeRandDbMat(1, n, 0, 0);
            double[,] net_ab = MakeRandDbMat(1, n, 0, 0);
            //%节点初始化
            //y=zeros(1,N);
            //net=zeros(1,n);
            //net_ab=zeros(1,n);
            double[,] d_Wjk = MakeRandDbMat(n, M, 0, 0);
            double[,] d_Wij = MakeRandDbMat(N, n, 0, 0);
            double[,] d_a = MakeRandDbMat(1, n, 0, 0);
            double[,] d_b = MakeRandDbMat(1, n, 0, 0);
            //%权值学习增量初始化
            //d_Wjk=zeros(n,M);
            //d_Wij=zeros(N,n);
            //d_a=zeros(1,n);
            //d_b=zeros(1,n);
            double[,] inputps = GetDbMatMaxMin(input);
            double[,] inputn = DbMapMinMax(input, inputps);
            double[] outputps = GetDbMatMaxMin(output);
            double[,] outputn = DbMapMinMax(output, outputps);
            //%% 输入输出数据归一化


            //[inputn,inputps]=mapminmax(input');
            //[outputn,outputps]=mapminmax(output'); 
            //inputn=inputn';
            //outputn=outputn';
            double[] error = new double[maxgen];
            for (int i = 0; i < maxgen; i++)
            {
                error[i] = 0;
                for (int kk = 0; kk < input.ElementAt(0).Count; kk++)
                {
                    double[] x = new double[input.Count];
                    double yqw = new double();
                    for (int j = 0; j < input.Count; j++)
                    {
                        x[j] = inputn[kk, j];
                    }
                    yqw = outputn[kk, 0];
                    for (int j = 0; j < n; j++)
                    {
                        for (int k = 0; k < M; k++)
                        {
                            net[0, j] = net[0, j] + Wjk[j, k] * x[k];
                            net_ab[0, j] = (net_ab[0, j] - b[0, j]) / a[0, j];
                        }
                        double temp = mymorlet(net_ab[0, j]);
                        for (int k = 0; k < N; k++)
                        {
                            y = y + Wij[k, j] * temp;
                        }
                    }
                    error[i] = error[i] + Math.Abs(yqw - y);

                    for (int j = 0; j < n; j++)
                    {
                        double temp = mymorlet(net_ab[0, j]);
                        for (int k = 0; k < N; k++)
                        {
                            d_Wij[k, j] = d_Wij[k, j] - (yqw - y) * temp;
                        }
                        temp = d_mylorlet(net_ab[0, j]);
                        for (int k = 0; k < M; k++)
                        {
                            for (int l = 0; l < N; l++)
                            {
                                d_Wjk[j, k] = d_Wjk[j, k] + (yqw - y) * Wij[l, j];
                            }
                            d_Wjk[j, k] = -d_Wjk[j, k] * temp * x[k] / a[0, j];
                        }
                        d_b[0, j] = d_b[0, j] + (yqw - y) * Wij[0, j];
                        d_b[0, j] = d_b[0, j] * temp / a[0, j];
                        d_a[0, j] = d_a[0, j] + (yqw - y) * Wij[0, j];
                        d_a[0, j] = d_a[0, j] * temp * ((net[0, j] - b[0, j]) / b[0, j]) / a[0, j];
                    }

                    //权值参数更新
                    for (int j = 0; j < N; j++)
                    {
                        for (int k = 0; k < n; k++)
                        {
                            Wij[j, k] = Wij[j, k] - lr1 * d_Wij[j, k];
                        }
                    }
                    for (int j = 0; j < n; j++)
                    {
                        for (int k = 0; k < M; k++)
                        {
                            Wjk[j, k] = Wjk[j, k] - lr1 * d_Wjk[j, k];
                        }
                    }
                    for (int j = 0; j < 1; j++)
                    {
                        for (int k = 0; k < n; k++)
                        {
                            b[j, k] = b[j, k] - lr2 * d_b[j, k];
                        }
                    }
                    for (int j = 0; j < 1; j++)
                    {
                        for (int k = 0; k < n; k++)
                        {
                            a[j, k] = a[j, k] - lr2 * d_a[j, k];
                        }
                    }
                    d_Wjk = MakeRandDbMat(n, M, 0, 0);
                    d_Wij = MakeRandDbMat(N, n, 0, 0);
                    d_a = MakeRandDbMat(1, n, 0, 0);
                    d_b = MakeRandDbMat(1, n, 0, 0);

                    y = 0;
                    net = MakeRandDbMat(1, n, 0, 0);
                    net_ab = MakeRandDbMat(1, n, 0, 0);

                    Wjk_1 = Wjk;
                    Wjk_2 = Wjk;
                    Wij_1 = Wij;
                    Wij_2 = Wij;
                    a_1 = a;
                    a_2 = a;
                    b_1 = b;
                    b_2 = b;
                }
            }
            double[,] xx = DbMapMinMax(input_test, inputps);
            double[] yuce = new double[input_test.ElementAt(0).Count];
            for (int i = 0; i < input_test.ElementAt(0).Count; i++)
            {
                double[] xx_test = new double[input_test.ElementAt(0).Count];
                for (int j = 0; j < input_test.Count; j++)
                {
                    xx_test[j] = input_test.ElementAt(j).ElementAt(i);
                }
                for (int j = 0; j < n; j++)
                {
                    for (int k = 0; k < M; k++)
                    {
                        net_ab[0, j] = net[0, j] + Wjk[j, k] * xx_test[k];
                        net_ab[0, j] = (net_ab[0, j] - b[0, j]) / a[0, j];
                    }
                    double temp = mymorlet(net_ab[0, j]);
                    for (int k = 0; k < N; k++)
                    {
                        y = y + Wij[k, j] * temp;
                    }
                }
                yuce[i] = y;
                y = 0.0;
                net = MakeRandDbMat(1, n, 0, 0);
                net_ab = MakeRandDbMat(1, n, 0, 0);
            }
            double[] yun = ReserveMaxMin(yuce, outputps);
            DataTable result = new DataTable();
            result.Columns.Add("时间", typeof(string));
            result.Columns.Add("数据", typeof(double));
            for (int i = 0; i < yun.Length; i++)
            {
                result.Rows.Add("预测",yun[i]);
            }
            return result.Copy();
        }

    }
}
