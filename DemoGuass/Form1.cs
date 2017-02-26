using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Windows.Forms;

namespace DemoGuass
{
    public partial class Form1 : Form
    {
        Log logtxt = null;
        Random rand_seed = null;
        private Queue<Queue<double>> input = null;//各个序列以列为Queue单位
        private Queue<Queue<double>> input_test = null;
        private Queue<double> output = null;
        private Queue<double> output_test = null;
        private Queue<double> source_data = null;
        private int items_count = new int();
        public int input_row = new int();
        public int input_test_row = new int();
        public int input_col = new int();
        public int sample_num = new int();//sample_row应该等于input_row
        public int forecast_num = new int();

        public bool IsFilled = false;
        public Form1()
        {
            InitializeComponent();
        }

        private void Form1_Load(object sender, EventArgs e)
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
            logtxt = new Log(AppDomain.CurrentDomain.BaseDirectory + @"/log/Log.txt");
            input = new Queue<Queue<double>>();
            input_test = new Queue<Queue<double>>();
            output = new Queue<double>();
            output_test = new Queue<double>();
            source_data = new Queue<double>();
            items_count = 0;
            input_row = 0;
            input_test_row = 0;
            input_col = 0;
            sample_num = 0;
            forecast_num = 0;
        }


        private void MakeRandMat(ref Queue<Queue<double>> mat, int row, int col, double min = 0.0, double max = 1.0)
        {
            Random rand = new Random(rand_seed.Next());
            int colDvalue = mat.Count - col;
            if (colDvalue >= 0)
            {
                for (int i = 0; i < colDvalue; i++)
                {
                    mat.Dequeue();
                }
                for (int i = 0; i < mat.Count; i++)
                {
                    int rowDvalue = mat.ElementAt(i).Count - row;
                    if (rowDvalue > 0)
                    {
                        for (int j = 0; j < rowDvalue; j++)
                        {
                            mat.ElementAt(i).Dequeue();
                        }
                    }
                    else
                    {
                        for (int j = 0; j < -rowDvalue; j++)
                        {
                            mat.ElementAt(i).Enqueue(rand.Next((int)min * 1000, (int)max * 1000) / 1000.0);
                        }
                    }

                }
            }
            else
            {
                for (int i = 0; i < -colDvalue; i++)
                {
                    Queue<double> tempQueue = new Queue<double>();
                    for (int j = 0; j < row; j++)
                    {
                        tempQueue.Enqueue(rand.Next((int)min * 1000, (int)max * 1000) / 1000.0);
                    }
                    mat.Enqueue(tempQueue);
                }
            }
        }

        private Queue<Queue<double>> MakeRandMat(int row, int col, double min = 0.0, double max = 1.0)
        {
            Random rand = new Random(rand_seed.Next());
            Queue<Queue<double>> mat = new Queue<Queue<double>>();
            for (int i = 0; i < col; i++)
            {
                Queue<double> tempQueue = new Queue<double>();
                for (int j = 0; j < row; j++)
                {
                    tempQueue.Enqueue(rand.Next((int)min * 1000, (int)max * 1000) / 1000.0);
                }
                mat.Enqueue(tempQueue);
            }
            return mat;
        }

        private double[,] MakeRandDbMat(int row, int col, double min = 0.0, double max = 1.0)
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

        private double[,] GetDbMatMaxMin(Queue<Queue<double>> mat)//获取矩阵列最大值最小值，最大在前，最小在后
        {
            double[,] result = new double[2, mat.Count];
            for (int i = 0; i < mat.Count; i++)
            {
                result[0, i] = mat.ElementAt(i).Max();
                result[1, i] = mat.ElementAt(i).Min();
            }
            return result;
        }

        private double[] GetDbMatMaxMin(Queue<double> mat)//获取矩阵列最大值最小值，最大在前，最小在后
        {
            double[] result = new double[2];
            result[0] = mat.Max();
            result[1] = mat.Min();
            return result;
        }

        private double[,] DbMapMinMax(Queue<double> mat, double[] maxmin)//二维矩阵归一化,返回数组
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

        private double[,] DbMapMinMax(Queue<Queue<double>> mat, double[,] maxmin)//二维矩阵归一化,返回数组
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

        private double[,] ReserveMaxMin(double[,] mat, double[,] maxmin)
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

        private double[] ReserveMaxMin(double[] mat, double[] maxmin)
        {
            int row = mat.Length;
            double[] result = new double[row];
            for (int i = 0; i < row; i++)
            {
                result[i] = (mat[i] + 1) * (maxmin[0] - maxmin[1]) / 2 + maxmin[1];
            }
            return result;
        }
        //private double[,] ApplyMapMinMax(Queue<Queue<double>> mat,double[,] maxmin)
        //{
        //    double[,] result = new double[mat.ElementAt(0).Count, mat.Count];
        //    for (int i = 0; i < mat.Count; i++)
        //    {
        //        for (int j = 0; j < mat.ElementAt(0).Count; j++)
        //        {
        //            if (maxmin[0, i] == maxmin[1, i])
        //                result[j, i] = 0.0;
        //            else
        //                result[j, i] = ((mat.ElementAt(i).ElementAt(j) - maxmin[1, i]) * 2) / (maxmin[0, i] - maxmin[1, i]) - 1;
        //        }
        //    }
        //    return result;
        //}
        private Queue<Queue<double>> GetMatMaxMin(Queue<Queue<double>> mat)//获取矩阵列最大值最小值，最大在前，最小在后
        {
            Queue<Queue<double>> result = new Queue<Queue<double>>();
            for (int i = 0; i < mat.Count; i++)
            {
                Queue<double> tempqueue = new Queue<double>();
                tempqueue.Enqueue(mat.ElementAt(i).Max());
                tempqueue.Enqueue(mat.ElementAt(i).Min());
                result.Enqueue(tempqueue);
            }
            return result;
        }



        private Queue<double> GetMatMaxMin(Queue<double> mat)//获取矩阵列最大值最小值，最大在前，最小在后
        {
            Queue<double> result = new Queue<double>();
            result.Enqueue(mat.Max());
            result.Enqueue(mat.Min());
            return result;
        }

        private Queue<Queue<double>> MapMinMax(Queue<Queue<double>> mat, Queue<Queue<double>> maxmin)//二维矩阵归一化
        {
            Queue<Queue<double>> result = new Queue<Queue<double>>();
            for (int i = 0; i < mat.Count; i++)
            {
                Queue<double> tempqueue = new Queue<double>();
                for (int j = 0; j < mat.ElementAt(0).Count; j++)
                {
                    if (maxmin.ElementAt(i).ElementAt(0) == maxmin.ElementAt(i).ElementAt(1))
                        tempqueue.Enqueue(new double());
                    else
                        tempqueue.Enqueue(((mat.ElementAt(i).ElementAt(j) - maxmin.ElementAt(i).ElementAt(1)) * 2) / (maxmin.ElementAt(i).ElementAt(0) - maxmin.ElementAt(i).ElementAt(1)) - 1);
                }
                result.Enqueue(tempqueue);
            }
            return result;
        }

        private Queue<double> MapMinMax(Queue<double> mat, Queue<double> maxmin)//二维矩阵归一化
        {
            Queue<double> result = new Queue<double>();
            for (int i = 0; i < mat.Count; i++)
            {
                //Queue<double> tempqueue = new Queue<double>();
                if (maxmin.ElementAt(0) == maxmin.ElementAt(1))
                    result.Enqueue(new double());
                else
                    result.Enqueue(((mat.ElementAt(i) - maxmin.ElementAt(1)) * 2) / (maxmin.ElementAt(0) - maxmin.ElementAt(1)) - 1);

            }
            return result;
        }

        private bool AddDataToSource(double data)
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

        public bool MakeQueue()
        {
            //while (input.Count<input_col)
            //{
            //    input.Enqueue(new  Queue<double>());
            //}
            //while (input_test.Count < input_col)
            //{
            //    input_test.Enqueue(new Queue<double>());
            //}
            //for (int i = 0; i < input_col; i++)
            //{
            //    while (input.ElementAt(i).Count > input_row)
            //    {
            //        input.ElementAt(i).Dequeue();
            //    }
            //    while (input_test.ElementAt(i).Count > input_test_row)
            //    {
            //        input_test.ElementAt(i).Dequeue();
            //    }

            //    while (input.ElementAt(i).Count<input_row)
            //    {
            //        input.ElementAt(i).Enqueue(new double());
            //    }
            //    while (input_test.ElementAt(i).Count < input_test_row)
            //    {
            //        input_test.ElementAt(i).Enqueue(new double());
            //    }
            //}
            //int items = input.Count * input.ElementAt(0).Count + input_test.Count * input_test.ElementAt(0).Count;
            if (items_count >= (input_row + input_test_row) * input_col)
            {
                input.Clear();
                input_test.Clear();
                output.Clear();
                for (int i = 0; i < items_count; )
                {
                    for (int j = 0; j < input_col; j++)
                    {
                        input.Enqueue(new Queue<double>());
                        for (int k = 0; k < input_row; k++)
                        {
                            input.ElementAt(j).Enqueue(source_data.ElementAt(items_count - 1 - i++));
                        }
                        input_test.Enqueue(new Queue<double>());
                        for (int k = 0; k < input_test_row; k++)
                        {
                            input_test.ElementAt(j).Enqueue(source_data.ElementAt(items_count - 1 - i++));
                        }
                    }
                }
                return true;
            }
            else
                return false;
        }

        private void AddDataToQueue(double data)
        {
            while (input.Count < input_col)
            {
                input.Enqueue(new Queue<double>());
            }
            while (input_test.Count < input_col)
            {
                input_test.Enqueue(new Queue<double>());
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
                    for (int j = 0; j < input_col; j++)
                    {
                        input.Enqueue(new Queue<double>());
                        for (int k = 0; k < input_row; k++)
                        {
                            input.ElementAt(j).Enqueue(source_data.ElementAt(FirstPos + i++));
                        }
                        input_test.Enqueue(new Queue<double>());
                        for (int k = 0; k < input_test_row; k++)
                        {
                            input_test.ElementAt(j).Enqueue(source_data.ElementAt(FirstPos + i++));
                        }
                    }



                    //int i = 0;
                    //for (int j = 0; j < input_col; j++)
                    //{
                    //    input.Enqueue(new Queue<double>());
                    //    for (int k = 0; k < input_row; k++)
                    //    {
                    //        input.ElementAt(j).Enqueue(source_data.ElementAt(items_count - 1 - i++));
                    //    }
                    //    input_test.Enqueue(new Queue<double>());
                    //    for (int k = 0; k < input_test_row; k++)
                    //    {
                    //        input_test.ElementAt(j).Enqueue(source_data.ElementAt(items_count - 1 - i++));
                    //    }
                    //}
                }
                else
                {
                    input.Clear();
                    input_test.Clear();
                    output.Clear();
                    int i = 0, l = 0;
                    int ZeroCount = (input_row + input_test_row) * input_col - items_count;
                    for (int j = 0; j < input_col; j++)
                    {
                        input.Enqueue(new Queue<double>());
                        for (int k = 0; k < input_row; k++)
                        {
                            if (l < ZeroCount)
                            {
                                input.ElementAt(j).Enqueue(new double());
                                l++;
                            }
                            else
                                input.ElementAt(j).Enqueue(source_data.ElementAt(i++));
                        }
                        input_test.Enqueue(new Queue<double>());
                        for (int k = 0; k < input_test_row; k++)
                        {
                            if (l < ZeroCount)
                            {
                                input_test.ElementAt(j).Enqueue(new double());
                                l++;
                            }
                            else
                                input_test.ElementAt(j).Enqueue(source_data.ElementAt(i++));
                        }
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
            //else if (input_row < input.ElementAt(0).Count || input_test_row < input_test.ElementAt(0).Count)
            //{
            //    input.Clear();
            //    input_test.Clear();
            //    output.Clear();
            //    int i = 0;
            //    if (items_count > (input_row + input_test_row) * input_col)
            //    {
            //        int FirstPos = items_count - (input_row + input_test_row) * input_col;
            //        for (int j = 0; j < input_col; j++)
            //        {
            //            input.Enqueue(new Queue<double>());
            //            for (int k = 0; k < input_row; k++)
            //            {
            //                input.ElementAt(j).Enqueue(source_data.ElementAt(FirstPos + i++));
            //            }
            //            input_test.Enqueue(new Queue<double>());
            //            for (int k = 0; k < input_test_row; k++)
            //            {
            //                input_test.ElementAt(j).Enqueue(source_data.ElementAt(FirstPos + i++));
            //            }
            //        }
            //        output.Clear();
            //        for (i = 0; i < sample_num; i++)
            //        {
            //            if (items_count >= sample_num)
            //                output.Enqueue(source_data.ElementAt(items_count - sample_num + i));
            //            else
            //            {
            //                if (i < sample_num - items_count)
            //                {
            //                    output.Enqueue(new double());
            //                }
            //                else
            //                    output.Enqueue(source_data.ElementAt(items_count - sample_num + i));
            //            }
            //        }

            //    }
            //    else
            //    {

            //    }
            //}
            //if (input.ElementAt(0).Count < input_row ||input_test.ElementAt(0).Count < input_test_row)
            //{

            //}
            //for (int i = 0; i < input_col; i++)
            //{
            //    while (input.ElementAt(i).Count > input_row)
            //    {
            //        input.ElementAt(i).Dequeue();
            //    }
            //    while (input_test.ElementAt(i).Count > input_test_row)
            //    {
            //        input_test.ElementAt(i).Dequeue();
            //    }

            //    while (input.ElementAt(i).Count < input_row)
            //    {
            //        input.ElementAt(i).Enqueue(new double());
            //    }
            //    while (input_test.ElementAt(i).Count < input_test_row)
            //    {
            //        input_test.ElementAt(i).Enqueue(new double());
            //    }
            //}
            //while (output.Count<sample_num)
            //{
            //    output.Enqueue(new double());  
            //}

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
        //public void MakeInputTestQueue(int row, int col)
        //{
        //}

        //public void ResizeOutputQueue(int length)
        //{

        //}

        private double mymorlet(double data)
        {
            return Math.Exp(-Math.Pow(data, 2) / 2) * Math.Cos(1.75 * data);
        }

        private double d_mylorlet(double data)
        {
            return -1.75 * Math.Sin(1.75 * data) * Math.Exp(-Math.Pow(data, 2) / 2) - data * Math.Cos(1.75 * data) * Math.Exp(-Math.Pow(data, 2) / 2);
        }

        void Forecast()
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

            //for (int i = 0; i < maxgen; i++)
            //{
            //    Queue<double> error = new Queue<double>();
            //    error.Enqueue(new double());
            //    for (int kk = 0; kk < input.ElementAt(0).Count; kk++)
            //    {
            //        Queue<double> x = new Queue<double>();
            //        Queue<double> yqw = new Queue<double>();
            //        for (int j = 0; j < inputn.Count; j++)
            //        {
            //            x.Enqueue(inputn.ElementAt(j).ElementAt(kk));
            //        }
            //        yqw.Enqueue(outputn.ElementAt(kk));
            //        Queue<Queue<double>> tempnet = new Queue<Queue<double>>();
            //        Queue<Queue<double>> tempnet_ab = new Queue<Queue<double>>();
            //        for (int j = 0; j < n; j++)
            //        {
            //            tempnet.Enqueue(new Queue<double>());
            //            tempnet_ab.Enqueue(new Queue<double>());
            //            double tempnetval = 0.0;
            //            double tempnet_abval = 0.0;
            //            for (int k = 0; k < M; k++)
            //            {
            //                tempnetval = net.ElementAt(j).ElementAt(0) + Wjk.ElementAt(k).ElementAt(j) * x.ElementAt(k);
            //                tempnet_abval = tempnetval - b.ElementAt(j).ElementAt(0) / a.ElementAt(j).ElementAt(0);
            //            }
            //            tempnet.ElementAt(j).Enqueue(tempnetval);
            //            tempnet_ab.ElementAt(j).Enqueue(tempnet_abval);
            //            double temp = mymorlet(tempnet_abval);
            //            Queue<Queue<double>> tempy = new Queue<Queue<double>>();
            //            for (int k = 0; k < N; k++)
            //            {
            //                tempy.Enqueue(new Queue<double>());
            //                tempy.ElementAt(k).Enqueue(y.ElementAt(k).ElementAt(0)+Wjk.ElementAt(j).ElementAt(k)*temp);
            //            }
            //            y = tempy;
            //        }
            //        net = tempnet;
            //        net_ab = tempnet_ab;
            //        double temperror = error.ElementAt(i)+;
            //        for (int j = 0; j < length; j++)
            //{

            //}
            //        temperror+=
            //    }
            //}


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
                        d_b[0, j] = d_b[0, j] + (yqw - y) * Wij[kk, j];
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
        }
            // 网络预测
            //预测输入归一化

            //%% 网络训练
            //for i=1:maxgen

            //    %误差累计
            //    error(i)=0;

            //    % 循环训练
            //    for kk=1:size(input,1)
            //        x=inputn(kk,:);
            //        yqw=outputn(kk,:);

            //        for j=1:n
            //            for k=1:M
            //                net(j)=net(j)+Wjk(j,k)*x(k);
            //                net_ab(j)=(net(j)-b(j))/a(j);
            //            end
            //            temp=mymorlet(net_ab(j));
            //            for k=1:N
            //                y=y+Wij(k,j)*temp;   %小波函数
            //            end
            //        end

            //        %计算误差和
            //        error(i)=error(i)+sum(abs(yqw-y));

            //        %权值调整
            //        for j=1:n
            //            %计算d_Wij
            //            temp=mymorlet(net_ab(j));
            //            for k=1:N
            //                d_Wij(k,j)=d_Wij(k,j)-(yqw(k)-y(k))*temp;
            //            end
            //            %计算d_Wjk
            //            temp=d_mymorlet(net_ab(j));
            //            for k=1:M
            //                for l=1:N
            //                    d_Wjk(j,k)=d_Wjk(j,k)+(yqw(l)-y(l))*Wij(l,j) ;
            //                end
            //                d_Wjk(j,k)=-d_Wjk(j,k)*temp*x(k)/a(j);
            //            end
            //            %计算d_b
            //            for k=1:N
            //                d_b(j)=d_b(j)+(yqw(k)-y(k))*Wij(k,j);
            //            end
            //            d_b(j)=d_b(j)*temp/a(j);
            //            %计算d_a
            //            for k=1:N
            //                d_a(j)=d_a(j)+(yqw(k)-y(k))*Wij(k,j);
            //            end
            //            d_a(j)=d_a(j)*temp*((net(j)-b(j))/b(j))/a(j);
            //        end

            //        %权值参数更新      
            //        Wij=Wij-lr1*d_Wij;
            //        Wjk=Wjk-lr1*d_Wjk;
            //        b=b-lr2*d_b;
            //        a=a-lr2*d_a;

            //        d_Wjk=zeros(n,M);
            //        d_Wij=zeros(N,n);
            //        d_a=zeros(1,n);
            //        d_b=zeros(1,n);

            //        y=zeros(1,N);
            //        net=zeros(1,n);
            //        net_ab=zeros(1,n);

            //        Wjk_1=Wjk;Wjk_2=Wjk_1;
            //        Wij_1=Wij;Wij_2=Wij_1;
            //        a_1=a;a_2=a_1;
            //        b_1=b;b_2=b_1;
            //    end
            //end

            //%% 网络预测
            //%预测输入归一化
            //x=mapminmax('apply',input_test',inputps);
            //x=x';

            //%网络预测
            //for i=1:92
            //    x_test=x(i,:);

            //    for j=1:1:n
            //        for k=1:1:M
            //            net(j)=net(j)+Wjk(j,k)*x_test(k);
            //            net_ab(j)=(net(j)-b(j))/a(j);
            //        end
            //        temp=mymorlet(net_ab(j));
            //        for k=1:N
            //            y(k)=y(k)+Wij(k,j)*temp ; 
            //        end
            //    end

            //    yuce(i)=y(k);
            //    y=zeros(1,N);
            //    net=zeros(1,n);
            //    net_ab=zeros(1,n);
            //end
            //%预测输出反归一化
            //ynn=mapminmax('reverse',yuce,outputps);

            //%% 结果分析
            //figure(1)
            //plot(ynn,'r*:')
            //hold on
            //plot(output_test,'bo--')
            //title('预测交通流量','fontsize',12)
            //legend('预测交通流量','实际交通流量')
            //xlabel('时间点')
            //ylabel('交通流量')

            //% web browser


        //            List<List<double>> Input = new List<List<double>>();
        //            List<List<double>> Output = new List<List<double>>();
        //            List<List<double>> PreInput = new List<List<double>>();
        //            List<List<double>> PreOutput = new List<List<double>>();
        //            int input_col;//定义训练输入样本列，数据计算以列为单位，与MATLAB中列向量同意义
        //            int input_row;//定义训练输入样本行
        //            int inputtest_row;//定义训练输出样本行
        //            int output_col;//定义预测输入样本列，即需要预测的数据，1列
        //            int output_row;//定义预测输入样本行
        //            int outputtest_row;//定义预测输出样本行，1列
        //            int maxgen; //迭代次数
        //            int n; //隐形节点个数	
        //            //double[,] Input = new double[20,1];
        //            ///////////////以下为个训练向量与预测向量结构初始化，预测时需进行数据有效化填充//////////////////////////////// 
        //            int i = 0;
        //            //Input = new List<double>(input_row);
        //            //Input = new List<List<double>>;
        //            Input = new List<List<double>>(input_row);
        //            for (i = 0; i < input_row; i++)
        //            {
        //                Input[i] = new List<double>(input_col);
        //            }
        //            Output = new List<List<double>>(inputtest_row);
        //            for (i = 0; i < inputtest_row; i++)
        //            {
        //                Output[i] = new List<double>(input_col);
        //            }
        //            PreInput = new List<List<double>>(output_row);
        //            for (i = 0; i < output_row; i++)
        //            {
        //                PreInput[i] = new List<double>(output_col);
        //            }
        //            PreOutput = new List<List<double>>(outputtest_row);
        //            for (i = 0; i < outputtest_row; i++)
        //            {
        //                PreOutput[i] = new List<double>(output_col);
        //            }



        //            int M = input_col;//输入节点个数
        //            int N = output_col;//输出节点个数
        //            double lr1 = 0.01; ////学习概率
        //            double lr2 = 0.001; ////学习概率
        //            int j, k, kk, kkk;
        //            //srand((unsigned)time( NULL ));
        //            double value = 0.0;
        //            string str = "";
        //            List<List<double>> Wjk = new List<List<double>>(n);
        //            for (i = 0; i < n; i++)
        //            {
        //                Wjk[i] = new List<double>(input_col);
        //            }
        //            List<List<double>> Wjk_1 = new List<List<double>>(n);
        //            for (i = 0; i < n; i++)
        //            {
        //                Wjk_1[i] = new List<double>(input_col);
        //            }
        //            List<List<double>> Wjk_2 = new List<List<double>>(n);
        //            for (i = 0; i < n; i++)
        //            {
        //                Wjk_2[i] = new List<double>(input_col);
        //            }

        //            for (i = 0; i < n; i++)
        //            {
        //                for (j = 0; j < input_col; j++)
        //                {
        //                    Wjk[i][j] = randn(2);
        //                    Wjk_1[i][j] = Wjk[i][j];
        //                    Wjk_2[i][j] = Wjk_1[i][j];
        //                }
        //            }
        //            List<List<double>> Wij = new List<List<double>>(output_col);
        //            for (i = 0; i < output_col; i++)
        //            {
        //                Wij[i] = new List<double>(n);
        //            }
        //            List<List<double>> Wij_1 = new List<List<double>>(output_col);
        //            for (i = 0; i < output_col; i++)
        //            {
        //                Wij_1[i] = new List<double>(n);
        //            }
        //            List<List<double>> Wij_2 = new List<List<double>>(output_col);
        //            for (i = 0; i < output_col; i++)
        //            {
        //                Wij_2[i] = new List<double>(n);
        //            }


        //            for (i = 0; i < output_col; i++)
        //            {
        //                for (j = 0; j < n; j++)
        //                {
        //                    Wij[i][j] = randn(2);
        //                    Wij_1[i][j] = Wij[i][j];
        //                    Wij_2[i][j] = Wij_1[i][j];
        //                }
        //            }

        //            List<List<double>> a = new List<List<double>>(1);
        //            for (i = 0; i < 1; i++)
        //            {
        //                a[i] = new List<double>(n);
        //            }
        //            List<List<double>> a_1 = new List<List<double>>(1);
        //            for (i = 0; i < 1; i++)
        //            {
        //                a_1[i] = new List<double>(n);
        //            }
        //            List<List<double>> a_2 = new List<List<double>>(1);
        //            for (i = 0; i < 1; i++)
        //            {
        //                a_2[i] = new List<double>(n);
        //            }

        //            for (i = 0; i < 1; i++)
        //            {
        //                for (j = 0; j < n; j++)
        //                {
        //                    a[i][j] = randn(2);
        //                    a_1[i][j] = a[i][j];
        //                    a_2[i][j] = a_1[i][j];
        //                }
        //            }

        //            List<List<double>> b = new List<List<double>>(1);
        //            for (i = 0; i < 1; i++)
        //            {
        //                b[i] = new List<double>(n);
        //            }
        //            List<List<double>> b_1 = new List<List<double>>(1);
        //            for (i = 0; i < 1; i++)
        //            {
        //                b_1[i] = new List<double>(n);
        //            }
        //            List<List<double>> b_2 = new List<List<double>>(1);
        //            for (i = 0; i < 1; i++)
        //            {
        //                b_2[i] = new List<double>(n);
        //            }
        //            for (i = 0; i < 1; i++)
        //            {
        //                for (j = 0; j < n; j++)
        //                {
        //                    b[i][j] = randn(2);
        //                    b_1[i][j] = b[i][j];
        //                    b_2[i][j] = b_1[i][j];
        //                }
        //            }
        //            List<List<double>> y = new List<List<double>>(1);
        //            for (i = 0; i < 1; i++)
        //            {
        //                y[i] = new List<double>(output_col);
        //            }
        //            List<List<double>> net = new List<List<double>>(1);
        //            for (i = 0; i < 1; i++)
        //            {
        //                net[i] = new List<double>(n);
        //            }
        //            List<List<double>> net_ab = new List<List<double>>(1);
        //            for (i = 0; i < 1; i++)
        //            {
        //                net_ab[i] = new List<double>(n);
        //            }
        //            List<List<double>> d_Wjk = new List<List<double>>(n);
        //            for (i = 0; i < n; i++)
        //            {
        //                d_Wjk[i] = new List<double>(input_col);
        //            }
        //            List<List<double>> d_Wij = new List<List<double>>(output_col);
        //            for (i = 0; i < output_col; i++)
        //            {
        //                d_Wij[i] = new List<double>(n);
        //            }
        //            List<List<double>> d_a = new List<List<double>>(1);
        //            for (i = 0; i < 1; i++)
        //            {
        //                d_a[i] = new List<double>(n);
        //            }
        //            List<List<double>> d_b = new List<List<double>>(1);
        //            for (i = 0; i < 1; i++)
        //            {
        //                d_b[i] = new List<double>(n);
        //            }
        //            ////// 输入输出数据归一化
        //            double tempmax = 0.0, tempmin = 0.0;
        //            List<double> input_max = new List<double>(input_col);
        //            List<double> input_min = new List<double>(input_col);
        //            List<List<double>> Input1 = new List<List<double>>(input_row);
        //            for (i = 0; i < input_row; i++)
        //            {
        //                Input1[i] = new List<double>(input_col);
        //            }
        //            for (i = 0; i < input_col; i++)
        //            {
        //                tempmax = Input[0][i];
        //                tempmin = Input[0][i];
        //                for (j = 0; j < input_row; j++)
        //                {
        //                    if (Input[j][i] > tempmax)
        //                        tempmax = Input[j][i];
        //                    if (Input[j][i] < tempmin)
        //                        tempmin = Input[j][i];
        //                }
        //                input_max[i] = tempmax;
        //                input_min[i] = tempmin;
        //                for (j = 0; j < input_row; j++)
        //                {
        //                    if (input_max[i] == input_min[i])
        //                        Input1[j][i] = 0;
        //                    else
        //                        Input1[j][i] = 2 * (Input[j][i] - input_min[i]) / (input_max[i] - input_min[i]) - 1;
        //                }
        //            }
        //            List<double> output_max = new List<double>(output_col);
        //            List<double> output_min = new List<double>(output_col);
        //            List<List<double>> PreInput1 = new List<List<double>>(output_row);
        //            for (i = 0; i < output_row; i++)
        //            {
        //                PreInput1[i] = new List<double>(output_col);
        //            }
        //            for (i = 0; i < output_col; i++)
        //            {
        //                tempmax = PreInput[0][i];
        //                tempmin = PreInput[0][i];
        //                for (j = 0; j < output_row; j++)
        //                {
        //                    if (PreInput[j][i] > tempmax)
        //                        tempmax = PreInput[j][i];
        //                    if (PreInput[j][i] < tempmin)
        //                        tempmin = PreInput[j][i];
        //                }
        //                output_max[i] = tempmax;
        //                output_min[i] = tempmin;
        //                for (j = 0; j < output_row; j++)
        //                {
        //                    if (output_max[i] == output_min[i])
        //                        PreInput1[j][i] = 0;
        //                    else
        //                        PreInput1[j][i] = 2 * (PreInput[j][i] - output_min[i]) / (output_max[i] - output_min[i]) - 1;
        //                }
        //            }

        //            ////// 网络训练
        //            List<double> x = new List<double>(input_col);
        //            List<double> yqw = new List<double>(output_col);
        //            List<double> error = new List<double>(maxgen);
        //            double temp = 0.0;
        //            for (i = 0; i < maxgen; i++)
        //            {
        //                ////误差累计
        //                error[i] = 0.0;
        //                for (kk = 0; kk < input_row; kk++)
        //                {
        //                    for (kkk = 0; kkk < input_col; kkk++)
        //                    {
        //                        x[kkk] = Input1[kk][kkk];
        //                    }
        //                    for (kkk = 0; kkk < output_col; kkk++)
        //                    {
        //                        yqw[kkk] = PreInput1[kk][kkk];
        //                    }
        //                    for (j = 0; j < n; j++)
        //                    {
        //                        for (k = 0; k < input_col; k++)
        //                        {
        //                            net[0][j] = net[0][j] + Wjk[j][k] * x[k];
        //                            net_ab[0][j] = (net[0][j] - b[0][j]) / a[0][j];
        //                        }
        //                        temp = mymorlet(net_ab[0][j]);
        //                        for (k = 0; k < output_col; k++)
        //                        {
        //                            y[0][k] = y[0][k] + Wij[k][j] * temp;
        //                        }
        //                    }
        //                    for (j = 0; j < output_col; j++)
        //                    {
        //                        temp = temp + abs(yqw[j] - y[0][j]);
        //                    }
        //                    error[i] = error[i] + temp;
        //                    for (j = 0; j < n; j++)
        //                    {
        //                        temp = mymorlet(net_ab[0][j]);
        //                        for (k = 0; k < output_col; k++)
        //                        {
        //                            d_Wij[k][j] = d_Wij[k][j] - (yqw[k] - y[0][k]) * temp;
        //                        }
        //                        temp = d_mymorlet(net_ab[0][j]);
        //                        for (k = 0; k < input_col; k++)
        //                        {
        //                            for (kkk = 0; kkk < output_col; kkk++)
        //                            {
        //                                d_Wjk[j][k] = d_Wjk[j][k] + (yqw[kkk] - y[0][kkk]) * Wij[kkk][j];
        //                            }
        //                            d_Wjk[j][k] = -d_Wjk[j][k] * temp * x[k] / a[0][j];
        //                        }
        //                        for (k = 0; k < output_col; k++)
        //                        {
        //                            d_b[0][j] = d_b[0][j] + (yqw[k] - y[0][k]) * Wij[k][j];
        //                        }
        //                        d_b[0][j] = d_b[0][j] * temp / a[0][j];
        //                        for (k = 0; k < output_col; k++)
        //                        {
        //                            d_a[0][j] = d_a[0][j] + (yqw[k] - y[0][k]) * Wij[k][j];
        //                        }
        //                        d_a[0][j] = d_a[0][j] * temp * ((net[0][j] - b[0][j]) / b[0][j]) / a[0][j];
        //                    }

        //                    ///权值参数更新
        //                    for (j = 0; j < n; j++)
        //                    {
        //                        for (k = 0; k < input_col; k++)
        //                        {
        //                            Wjk[j][k] = Wjk[j][k] - lr1 * d_Wjk[j][k];
        //                            d_Wjk[j][k] = 0.0;
        //                            Wjk_1[j][k] = Wjk[j][k];
        //                            Wjk_2[j][k] = Wjk_1[j][k];
        //                        }
        //                        b[0][j] = b[0][j] - lr2 * d_b[0][j];
        //                        a[0][j] = a[0][j] - lr2 * d_a[0][j];
        //                        d_a[0][j] = 0.0;
        //                        d_b[0][j] = 0.0;
        //                        net[0][j] = 0.0;
        //                        net_ab[0][j] = 0.0;
        //                        a_1[0][j] = a[0][j];
        //                        a_2[0][j] = a_1[0][j];
        //                        b_1[0][j] = b[0][j];
        //                        b_1[0][j] = b_1[0][j];
        //                    }
        //                    for (j = 0; j < output_col; j++)
        //                    {
        //                        for (k = 0; k < n; k++)
        //                        {
        //                            Wij[j][k] = Wij[j][k] - lr1 * d_Wij[j][k];
        //                            d_Wij[j][k] = 0.0;
        //                            Wij_1[j][k] = Wij[j][k];
        //                            Wij_2[j][k] = Wij_1[j][k];
        //                        }
        //                        y[0][j] = 0.0;
        //                    }
        //                }
        //            }


        //            ////// 网络预测
        //            //////预测输入归一化
        //            List<double> inputtest_max = new List<double>(input_col);
        //            List<double> inputtest_min = new List<double>(input_col);
        //            List<List<double>> OutputAsInput1 = new List<List<double>>(inputtest_row);
        //            for (i = 0; i < inputtest_row; i++)
        //            {
        //                OutputAsInput1[i] = new List<double>(input_col);
        //            }
        //            //double inputtest_min[input_col] = {{0.0}};
        //            for (i = 0; i < input_col; i++)
        //            {
        //                tempmax = Output[0][i];
        //                tempmin = Output[0][i];
        //                for (j = 0; j < inputtest_row; j++)
        //                {
        //                    if (Output[j][i] > tempmax)
        //                        tempmax = Output[j][i];
        //                    if (Output[j][i] < tempmin)
        //                        tempmin = Output[j][i];
        //                }
        //                inputtest_max[i] = tempmax;
        //                inputtest_min[i] = tempmin;
        //                for (j = 0; j < inputtest_row; j++)
        //                {
        //                    if (inputtest_max[i] == inputtest_min[i])
        //                        OutputAsInput1[j][i] = 0.5 * (input_max[i] + input_min[i]);
        //                    else
        //                        OutputAsInput1[j][i] = (Output[j][i] - inputtest_min[i]) / (inputtest_max[i] - inputtest_min[i]) * (input_max[i] - input_min[i]) + input_min[i];
        //                }
        //            }

        //            /////网络预测
        //            List<double> x_test = new List<double>(input_col);
        //            List<List<double>> yuce = new List<List<double>>(1);
        //            for (i = 0; i < 1; i++)
        //            {
        //                yuce[i] = new List<double>(inputtest_row);
        //            }
        //            for (i = 0; i < inputtest_row; i++)
        //            {
        //                for (j = 0; j < input_col; j++)
        //                {
        //                    x_test[j] = OutputAsInput1[i][j];
        //                }
        //                for (j = 0; j < n; j++)
        //                {
        //                    for (k = 0; k < input_col; k++)
        //                    {
        //                        net[0][j] = net[0][j] + Wjk[j][k] * x_test[k];
        //                        net_ab[0][j] = (net[0][j] - b[0][j]) / a[0][j];
        //                    }
        //                    temp = mymorlet(net_ab[0][j]);
        //                    for (k = 0; k < output_col; k++)
        //                    {
        //                        y[0][k] = y[0][k] + Wij[k][j] * temp;
        //                    }
        //                }
        //                yuce[0][i] = y[0][k - 1];
        //                for (j = 0; j < output_col; j++)
        //                {
        //                    y[0][j] = 0.0;
        //                }
        //                for (j = 0; j < n; j++)
        //                {
        //                    net[0][j] = 0.0;
        //                    net_ab[0][j] = 0.0;
        //                }
        //            }

        //            List<double> yun = new List<double>(outputtest_row);
        //            for (i = 0; i < outputtest_row; i++)
        //            {
        //                yun[i] = (yuce[0][i] + 1) / 2 * (output_max[0] - output_min[0]) + output_min[0];
        //                PreOutput[i][output_col - 1] = yun[i];
        //            }
        //            return;
        //        }





        private void button1_Click(object sender, EventArgs e)
        {

            Queue<Queue<double>> mat = MakeRandMat(5, 7, 0, 0);
            for (int i = 0; i < mat.Count; i++)
            {
                string str1 = "";
                for (int j = 0; j < mat.ElementAt(i).Count; j++)
                {
                    str1 = str1 + "|" + mat.ElementAt(i).ElementAt(j).ToString();
                }
                logtxt.log(str1);
            }
            logtxt.log("0000000000000000000000000000000000000000000000000000000000000000");
            return;


            Random ran = new Random(rand_seed.Next());
            input_col = 4;
            input_row = 30;
            input_test_row = 10;
            sample_num = 20;
            AddDataToSource(ran.Next(0, 1000) / 1000.0);

            string str = "";
            for (int i = 0; i < input_row; i++)
            {

                str = string.Format("{0},{1},{2},{3}", input.ElementAt(0).ElementAt(i).ToString("f3"), input.ElementAt(1).ElementAt(i).ToString("f3"), input.ElementAt(2).ElementAt(i).ToString("f3"), input.ElementAt(3).ElementAt(i).ToString("f3"));
                logtxt.log(str);
            }
            logtxt.log("0000000000000000000000000000000000000000000000000000000000000000");
            for (int i = 0; i < input_test_row; i++)
            {

                str = string.Format("{0},{1},{2},{3}", input_test.ElementAt(0).ElementAt(i).ToString("f3"), input_test.ElementAt(1).ElementAt(i).ToString("f3"), input_test.ElementAt(2).ElementAt(i).ToString("f3"), input_test.ElementAt(3).ElementAt(i).ToString("f3"));
                logtxt.log(str);
            }
            logtxt.log("0000000000000000000000000000000000000000000000000000000000000000");
            for (int i = 0; i < sample_num; i++)
            {
                str = string.Format("{0}", output.ElementAt(i).ToString("f3"));
                logtxt.log(str);
            }
            logtxt.log("0000000000000000000000000000000000000000000000000000000000000000");
            //MakeRandMat(ref input, 20, 30);
            //for (int i = 0; i < input.Count; i++)
            //{
            //    string str = "";
            //    for (int j = 0; j < input.ElementAt(i).Count; j++)
            //    {
            //        str = str + "|" + input.ElementAt(i).ElementAt(j).ToString();
            //    }
            //    logtxt.log(str);
            //}

        }

        private void button2_Click(object sender, EventArgs e)
        {
            double[,] result = MakeRandDbMat(3, 5, 0, 1);
            for (int i = 0; i < 3; i++)
            {
                string str = "";
                for (int j = 0; j < 5; j++)
                {
                    str = string.Format("{0},{1},{2},{3},{4}", result[i, 0], result[i, 1], result[i, 2], result[i, 3], result[i, 4]);
                }
                logtxt.log(str);
            }
            //Random ran = new Random(rand_seed.Next());
            //input_col = 4;
            //input_row = 50;
            //input_test_row = 20;
            //sample_num = 30;
            //AddDataToSource(ran.Next(0, 1000) / 1000.0);

            //string str = "";
            //for (int i = 0; i < input_row; i++)
            //{

            //    str = string.Format("{0},{1},{2},{3}", input.ElementAt(0).ElementAt(i).ToString("f3"), input.ElementAt(1).ElementAt(i).ToString("f3"), input.ElementAt(2).ElementAt(i).ToString("f3"), input.ElementAt(3).ElementAt(i).ToString("f3"));
            //    logtxt.log(str);
            //}
            //logtxt.log("0000000000000000000000000000000000000000000000000000000000000000");
            //for (int i = 0; i < input_test_row; i++)
            //{

            //    str = string.Format("{0},{1},{2},{3}", input_test.ElementAt(0).ElementAt(i).ToString("f3"), input_test.ElementAt(1).ElementAt(i).ToString("f3"), input_test.ElementAt(2).ElementAt(i).ToString("f3"), input_test.ElementAt(3).ElementAt(i).ToString("f3"));
            //    logtxt.log(str);
            //}
            //logtxt.log("0000000000000000000000000000000000000000000000000000000000000000");
            //for (int i = 0; i < sample_num; i++)
            //{
            //    str = string.Format("{0}", output.ElementAt(i).ToString("f3"));
            //    logtxt.log(str);
            //}
            //logtxt.log("0000000000000000000000000000000000000000000000000000000000000000");
            //Queue<Queue<double>> mat = MakeRandMat( 10, 5);
            //double total = 0.0;
            //for (int i = 0; i < mat.Count; i++)
            //{
            //    string str = "";
            //    for (int j = 0; j < mat.ElementAt(i).Count; j++)
            //    {
            //        total = total + mat.ElementAt(i).ElementAt(j);
            //        str = str + "|" + mat.ElementAt(i).ElementAt(j).ToString();
            //    }
            //    logtxt.log(str);
            //}
            //MakeRandMat(ref input, 10, 5);
            //double total = 0.0;
            //for (int i = 0; i < input.Count; i++)
            //{
            //    string str = "";
            //    for (int j = 0; j < input.ElementAt(i).Count; j++)
            //    {
            //        total = total + input.ElementAt(i).ElementAt(j);
            //        str = str + "|" + input.ElementAt(i).ElementAt(j).ToString();
            //    }
            //    logtxt.log(str);
            //}
            //logtxt.log(String.Format("{0}0000000000000000000000000000000000000000000000000000000000000000",total));
        }
    }
}
