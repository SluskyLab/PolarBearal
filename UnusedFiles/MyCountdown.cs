using System;
using System.Collections.Generic;
using System.Text;
using System.Collections;
using System.IO;
using System.Xml;
using System.Data;
using System.Runtime.Serialization.Formatters.Binary;
using System.Runtime.Serialization;
using System.Threading;

namespace betaBarrelProgram
{
    public class MyCountdown
    {
        object _locker = new object();
        int _value;

        public MyCountdown() { }
        public MyCountdown(int initialCount) { _value = initialCount; }

        public void Signal() { AddCount(-1); }

        public void AddCount(int amount)
        {
            lock (_locker)
            {
                _value += amount;
                if (_value <= 0) Monitor.PulseAll(_locker);
            }
        }

        public void Wait()
        {
            lock (_locker)
                while (_value > 0)
                    Monitor.Wait(_locker);
        }
    }

}
