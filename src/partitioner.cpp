#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cassert>
#include <vector>
#include <cmath>
#include <map>
#include "cell.h"
#include "net.h"
#include "partitioner.h"
using namespace std;

void Partitioner::parseInput(fstream &inFile)
{
    string str;
    // Set balance factor
    inFile >> str;
    _bFactor = stod(str);

    // Set up whole circuit
    while (inFile >> str)
    {
        if (str == "NET")
        {
            string netName, cellName, tmpCellName = "";
            inFile >> netName;
            int netId = _netNum;
            _netArray.push_back(new Net(netName));
            _netName2Id[netName] = netId;
            while (inFile >> cellName)
            {
                if (cellName == ";")
                {
                    tmpCellName = "";
                    break;
                }
                else
                {
                    // a newly seen cell
                    if (_cellName2Id.count(cellName) == 0)
                    {
                        int cellId = _cellNum;
                        _cellArray.push_back(new Cell(cellName, 0, cellId));
                        _cellName2Id[cellName] = cellId;
                        _cellArray[cellId]->addNet(netId);
                        _cellArray[cellId]->incPinNum();
                        _netArray[netId]->addCell(cellId);
                        ++_cellNum;
                        tmpCellName = cellName;
                    }
                    // an existed cell
                    else
                    {
                        if (cellName != tmpCellName)
                        {
                            assert(_cellName2Id.count(cellName) == 1);
                            int cellId = _cellName2Id[cellName];
                            _cellArray[cellId]->addNet(netId);
                            _cellArray[cellId]->incPinNum();
                            _netArray[netId]->addCell(cellId);
                            tmpCellName = cellName;
                        }
                    }
                }
            }
            ++_netNum;
        }
    }
    return;
}

void Partitioner::init_part()
{
    for (int i = 0; i < _cellNum; i++)
    {
        if (_cellArray[i]->getPinNum() > _maxPinNum)
            _maxPinNum++;
        if (i < _cellNum / 2)
        {
            _cellArray[i]->setPart(0);
            _partSize[0]++;
            vector<int> list = _cellArray[i]->getNetList();
            for (int j = 0; j < list.size(); j++)
            {
                _netArray[list[j]]->incPartCount(0);
            }
        }
        else
        {
            _cellArray[i]->setPart(1);
            _partSize[1]++;
            vector<int> list = _cellArray[i]->getNetList();
            for (int j = 0; j < list.size(); j++)
            {
                _netArray[list[j]]->incPartCount(1);
            }
        }
    }
}
void Partitioner::set_cutSize()
{
    _cutSize = 0;
    for (int i = 0; i < _netNum; i++)
    {
        if (_netArray[i]->getPartCount(0) > 0 && _netArray[i]->getPartCount(1) > 0)
            _cutSize++;
    }
}
void Partitioner::add_cell(int i)
{
    // cout << i << endl;
    map<int, Node *>::iterator itr;
    Node *node = _cellArray[i]->getNode();
    node->setNext(nullptr);
    node->setPrev(nullptr);
    // get gain
    int gain = _cellArray[i]->getGain();
    // cell in A(0)
    if (!_cellArray[i]->getPart())
    {
        itr = _bList[0].find(gain);
        if (itr == _bList[0].end())
        {
            _bList[0][gain] = node;
        }
        else
        {
            (*itr).second->setPrev(node);
            node->setNext((*itr).second);
            _bList[0][gain] = node;
        }
    }
    // cell in B(1)
    else
    {
        itr = _bList[1].find(gain);
        if (itr == _bList[1].end())
        {
            _bList[1][gain] = node;
        }
        else
        {
            (*itr).second->setPrev(node);
            node->setNext((*itr).second);
            _bList[0][gain] = node;
        }
    }
}
void Partitioner::init_gain()
{
    // clear b_list
    _bList[0].clear();
    _bList[1].clear();
    // go through every cell
    for (int i = 0; i < _cellNum; i++)
    {
        _cellArray[i]->setGain(0);
        vector<int> netlist = _cellArray[i]->getNetList();
        // calculate gain
        if (!_cellArray[i]->getPart())
        {
            for (int j = 0; j < netlist.size(); j++)
            {
                if (_netArray[netlist[j]]->getPartCount(0) == 1)
                    _cellArray[i]->incGain();
                if (_netArray[netlist[j]]->getPartCount(1) == 0)
                    _cellArray[i]->decGain();
            }
        }
        else
        {
            for (int j = 0; j < netlist.size(); j++)
            {
                if (_netArray[netlist[j]]->getPartCount(1) == 1)
                    _cellArray[i]->incGain();
                if (_netArray[netlist[j]]->getPartCount(0) == 0)
                    _cellArray[i]->decGain();
            }
        }
        // add cell
        add_cell(i);
        std::cout << _cellArray[i]->getName() << " " << _cellArray[i]->getGain() << " " << _cellArray[i]->getPinNum() << " " << _cellArray[i]->getPart() << endl;
    }
}
void Partitioner::update_gain(int i)
{
    vector<int> list = _cellArray[i]->getNetList();
    // check all nets
    for (int j = 0; j < list.size(); j++)
    {
        cout << _netArray[list[j]]->getName() << " " << _netArray[list[j]]->getPartCount(0) << " " << _netArray[list[j]]->getPartCount(1) << endl;
        vector<int> cells = _netArray[list[j]]->getCellList();

        // cell[i] in A(0)
        if (!_cellArray[i]->getPart())
        {
            // T(n)=0
            if (_netArray[list[j]]->getPartCount(1) == 0)
            {
                cout << "T(n)=0" << endl;
                for (int k = 0; k < cells.size(); k++)
                {
                    _cellArray[cells[k]]->incGain();
                    cout << _cellArray[cells[k]]->getName() << "++" << endl;
                }
            }
            // T(n)=1
            else if (_netArray[list[j]]->getPartCount(1) == 1)
            {
                cout << "T(n)=1" << endl;
                for (int k = 0; k < cells.size(); k++)
                {
                    if (_cellArray[cells[k]]->getPart())
                    {
                        _cellArray[cells[k]]->decGain();
                        cout << _cellArray[cells[k]]->getName() << "--" << endl;
                    }
                }
            }
            // lock cell[i]
            _cellArray[i]->lock();
            _netArray[list[j]]->incPartCount(1);
            _netArray[list[j]]->decPartCount(0);
            // F(n)=0

            if (_netArray[list[j]]->getPartCount(0) == 0)
            {
                cout << "F(n)=0" << endl;
                for (int k = 0; k < cells.size(); k++)
                {
                    _cellArray[cells[k]]->decGain();
                    cout << _cellArray[cells[k]]->getName() << "--" << endl;
                }
            }
            // F(n)=1

            else if (_netArray[list[j]]->getPartCount(0) == 1)
            {
                cout << "F(n)=1" << endl;
                for (int k = 0; k < cells.size(); k++)
                {
                    if (!_cellArray[cells[k]]->getPart())
                    {
                        _cellArray[cells[k]]->incGain();
                        cout << _cellArray[cells[k]]->getName() << "++" << endl;
                    }
                }
            }
        }
        // cell[i] in B(1)
        else
        {
            // T(n)=0
            if (_netArray[list[i]]->getPartCount(0) == 0)
            {
                for (int j = 0; j < cells.size(); j++)
                {
                    _cellArray[cells[j]]->incGain();
                }
            }
            // T(n)=1
            else if (_netArray[list[i]]->getPartCount(0) == 1)
            {
                for (int j = 0; j < cells.size(); j++)
                {
                    if (!_cellArray[cells[j]]->getPart())
                        _cellArray[cells[j]]->decGain();
                }
            }
            // lock cell[i]
            _cellArray[i]->lock();
            _netArray[list[i]]->incPartCount(0);
            _netArray[list[i]]->decPartCount(1);
            // F(n)=0
            if (_netArray[list[i]]->getPartCount(1) == 0)
            {
                for (int j = 0; j < cells.size(); j++)
                {
                    _cellArray[cells[j]]->decGain();
                }
            }
            // F(n)=1
            else if (_netArray[list[i]]->getPartCount(1) == 1)
            {
                for (int j = 0; j < cells.size(); j++)
                {
                    if (_cellArray[cells[j]]->getPart())
                        _cellArray[cells[j]]->incGain();
                }
            }
        }
    }
    // print gain
    for (int i = 0; i < _cellArray.size(); i++)
    {
        std::cout << _cellArray[i]->getName() << " " << _cellArray[i]->getGain() << endl;
    }
}
void Partitioner::partition()
{
    // initial partition
    init_part();
    set_cutSize();
    init_gain();
    update_gain(1);
    // cout << _netArray[3]->getPartCount(0) << endl;
    // cout << _netArray[3]->getPartCount(1) << endl;
    //  std::cout<<_cellArray[0]->getNetList()[];
}

void Partitioner::printSummary() const
{
    cout << endl;
    cout << "==================== Summary ====================" << endl;
    cout << " Cutsize: " << _cutSize << endl;
    cout << " Total cell number: " << _cellNum << endl;
    cout << " Total net number:  " << _netNum << endl;
    cout << " Cell Number of partition A: " << _partSize[0] << endl;
    cout << " Cell Number of partition B: " << _partSize[1] << endl;
    cout << "=================================================" << endl;
    cout << endl;
    return;
}

void Partitioner::reportNet() const
{
    cout << "Number of nets: " << _netNum << endl;
    for (size_t i = 0, end_i = _netArray.size(); i < end_i; ++i)
    {
        cout << setw(8) << _netArray[i]->getName() << ": ";
        vector<int> cellList = _netArray[i]->getCellList();
        for (size_t j = 0, end_j = cellList.size(); j < end_j; ++j)
        {
            cout << setw(8) << _cellArray[cellList[j]]->getName() << " ";
        }
        cout << endl;
    }
    return;
}

void Partitioner::reportCell() const
{
    cout << "Number of cells: " << _cellNum << endl;
    for (size_t i = 0, end_i = _cellArray.size(); i < end_i; ++i)
    {
        cout << setw(8) << _cellArray[i]->getName() << ": ";
        vector<int> netList = _cellArray[i]->getNetList();
        for (size_t j = 0, end_j = netList.size(); j < end_j; ++j)
        {
            cout << setw(8) << _netArray[netList[j]]->getName() << " ";
        }
        cout << endl;
    }
    return;
}

void Partitioner::writeResult(fstream &outFile)
{
    stringstream buff;
    buff << _cutSize;
    outFile << "Cutsize = " << buff.str() << '\n';
    buff.str("");
    buff << _partSize[0];
    outFile << "G1 " << buff.str() << '\n';
    for (size_t i = 0, end = _cellArray.size(); i < end; ++i)
    {
        if (_cellArray[i]->getPart() == 0)
        {
            outFile << _cellArray[i]->getName() << " ";
        }
    }
    outFile << ";\n";
    buff.str("");
    buff << _partSize[1];
    outFile << "G2 " << buff.str() << '\n';
    for (size_t i = 0, end = _cellArray.size(); i < end; ++i)
    {
        if (_cellArray[i]->getPart() == 1)
        {
            outFile << _cellArray[i]->getName() << " ";
        }
    }
    outFile << ";\n";
    return;
}

void Partitioner::clear()
{
    for (size_t i = 0, end = _cellArray.size(); i < end; ++i)
    {
        delete _cellArray[i];
    }
    for (size_t i = 0, end = _netArray.size(); i < end; ++i)
    {
        delete _netArray[i];
    }
    return;
}
