﻿using System;
using System.Collections.Generic;
using System.Text;
using System.IO;

namespace ClassForecast
{
    public class Log
    {
        private string logFile;
        private StreamWriter writer;
        private FileStream fileStream = null;

        public Log(string fileName)
        {
            logFile = fileName;
            CreateDirectory(logFile);
        }
        //使用
        //Log log = new Log(AppDomain.CurrentDomain.BaseDirectory + @"/log/Log.txt");
        //log.log(basePath);
        public void log(string info)
        {
            try
            {
                System.IO.FileInfo fileInfo = new System.IO.FileInfo(logFile);
                if (!fileInfo.Exists)
                {
                    fileStream = fileInfo.Create();
                    writer = new StreamWriter(fileStream);
                }
                else
                {
                    fileStream = fileInfo.Open(FileMode.Append, FileAccess.Write);
                    writer = new StreamWriter(fileStream);
                }
                writer.WriteLine(DateTime.Now.TimeOfDay.ToString() + ": " + info);
                //writer.WriteLine("----------------------------------");
            }
            finally
            {
                if (writer != null)
                {
                    writer.Close();
                    writer.Dispose();
                    fileStream.Close();
                    fileStream.Dispose();
                }
            }
        }

        public void CreateDirectory(string infoPath)
        {
            DirectoryInfo directoryInfo = Directory.GetParent(infoPath);
            if (!directoryInfo.Exists)
            {
                directoryInfo.Create();
            }
        }
    }
}