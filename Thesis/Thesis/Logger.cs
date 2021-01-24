using System;
using System.Collections.Generic;
using System.IO;
using System.Data;
using System.Text;
using System.Linq;

namespace Thesis
{
    class Logger : IDisposable
    {
        readonly StreamWriter writer;
        public readonly string filePath;

        /// <summary> Wraps a streamwriter for outputting content to a text file </summary>
        /// <remarks> Be sure to use a 'using' block or close the writer when it's no longer needed. </remarks>
        /// <param name="fname"> The name of the file you want it to create. A number will be appended if this name is already being used. </param>
        public Logger(string fname)
        {
            try
            {
                filePath = CreateUniqueOutputFilePath(fname);
                writer = new StreamWriter(filePath) { AutoFlush = true };
            }
            catch (Exception e) { Console.WriteLine(e); }
        }

        public void Write(string str) { writer?.Write(str); }
        public void WriteLine(string str) { writer?.WriteLine(str); }

        /// <summary> Writes the contents of a table to the log </summary>
        /// <remarks> Based on code found here https://stackoverflow.com/questions/4959722/c-sharp-datatable-to-csv  </remarks>
        public void WriteTable(DataTable table)
        {
            // --- Header ---
            WriteLine(table.TableName);

            // --- Contents ---
            StringBuilder builder = new StringBuilder();
            IEnumerable<string> columnNames = table.Columns.Cast<DataColumn>().
                                              Select(column => column.ColumnName);
            builder.AppendLine(string.Join(",", columnNames));

            foreach (DataRow row in table.Rows)
            {
                IEnumerable<string> fields = row.ItemArray.Select(field => field.ToString());
                builder.AppendLine(string.Join(",", fields));
            }

            writer.Write(builder.ToString());
        }

        public void WriteTablesSideBySide(params DataTable[] dataTables)
        {
            StringBuilder builder = new StringBuilder();
            // --- Header ---
            for (int i = 0; i < dataTables.Length; i++)
            {
                builder.Append(dataTables[i].TableName);
                builder.Append(',');
                for (int j = 0; j < dataTables[i].Columns.Count; j++)
                {
                    if (i < dataTables.Length - 1 || j < dataTables[i].Columns.Count - 1)
                        builder.Append(',');
                }
            }
            builder.AppendLine();


            // Write column names
            for (int tableIdx = 0; tableIdx < dataTables.Length; tableIdx++)
            {
                for (int col = 0; col < dataTables[tableIdx].Columns.Count; col++)
                {
                    builder.Append(dataTables[tableIdx].Columns[col].ColumnName + ",");
                }
                builder.Append(",");
            }
            builder.AppendLine();

            // Figure out how many rows we need
            int numRows = 0;
            for (int i = 0; i < dataTables.Length; i++)
            {
                numRows = Math.Max(numRows, dataTables[i].Rows.Count);
            }

            // --- Contents ---
            for (int rowIdx = 0; rowIdx < numRows; rowIdx++)
            {
                for (int tableIdx = 0; tableIdx < dataTables.Length; tableIdx++)
                {
                    // Cell Contents
                    if (rowIdx < dataTables[tableIdx].Rows.Count)
                    {
                        var row = dataTables[tableIdx].Rows[rowIdx];
                        IEnumerable<string> fields = row.ItemArray.Select(field => field.ToString());
                        builder.Append(string.Join(",", fields));
                        builder.Append(',');
                    }
                    else
                    {
                        // Fill in the space with empty cells
                        for (int col = 0; col < dataTables[tableIdx].Columns.Count; col++)
                        {
                            builder.Append(',');
                        }
                    }
                    builder.Append(','); // Add space between tables
                }
                if (rowIdx < numRows - 1) builder.AppendLine();
            }

            writer.Write(builder.ToString());
        }

        public static string CreateUniqueOutputFilePath(string fname)
        {
            // --- Append a unique number to the desired file name ---
            string currentDir = Environment.CurrentDirectory;
            string[] files = Directory.GetFiles(currentDir);
            List<string> names = new List<string>(files.Length);
            for (int i = 0; i < files.Length; i++)
            {
                names.Add(Path.GetFileName(files[i]));
            }
            int numberToAppend = 2;
            string fileName = fname;
            while (names.Contains(fileName))
            {
                fileName = fname.Insert(fname.LastIndexOf('.'), " " + numberToAppend.ToString());
                numberToAppend++;
            }
            
            return Environment.CurrentDirectory + Path.DirectorySeparatorChar + fileName;
        }


        public void Dispose()
        {
            writer?.Close();
        }
    }
}
