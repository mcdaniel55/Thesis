using System;
using System.Collections.Generic;
using System.Threading;
using System.Threading.Tasks;
using System.Diagnostics;

namespace ThesisOptNumericalTest
{
    static class ThreadTesting
    {
        public static void Test1()
        {
            Stopwatch watch = new Stopwatch();

            int number = 0;
            int size = 300000;
            int numTasks = 100;
            Action count = () => { for (int i = 0; i < size; i++) { if (i == size-1) Interlocked.Increment(ref number); } };
            Task[] tasks;

            for (int iteration = 0; iteration < 50; iteration++)
            {
                tasks = new Task[numTasks];
                for (int i = 0; i < tasks.Length; i++)
                {
                    tasks[i] = new Task(count);
                }
                watch.Restart();
                for (int i = 0; i < tasks.Length; i++)
                {
                    tasks[i].Start();
                }
                Task.WaitAll(tasks);
                watch.Stop();

                Console.WriteLine($"number = {number}");
                Console.WriteLine($"Time = {watch.Elapsed.TotalMilliseconds}ms");
            }

            number = 0;
            size = size * numTasks / 4;

            for (int iteration = 0; iteration < 50; iteration++)
            {
                tasks = new Task[4];
                
                for (int i = 0; i < tasks.Length; i++) { tasks[i] = new Task(count); }

                watch.Restart();
                for (int i = 0; i < tasks.Length; i++)
                {
                    tasks[i].Start();
                }
                Task.WaitAll(tasks);
                watch.Stop();

                Console.WriteLine($"number = {number}");
                Console.WriteLine($"Time = {watch.Elapsed.TotalMilliseconds}ms");
            }
            Console.ReadKey();
        }

        

    }
}
