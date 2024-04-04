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
    max = (1+_bFactor)*_cellNum/2;
    min = (1-_bFactor)*_cellNum/2;
    for (int i = 0; i < _cellNum; i++)
    {
        if (_cellArray[i]->getPinNum() > _maxPinNum)
            _maxPinNum++;
        if (i < _cellNum / 2)
        {
            _cellArray[i]->setPart(0);
        }
        else
        {
            _cellArray[i]->setPart(1);
        }
    }
}
void Partitioner::init_size_and_count(){
    //reset size and count
    _partSize[0] = 0;
    _partSize[1] = 0;
    for(int i = 0; i < _netNum; i++){
        _netArray[i]->setPartCount(0,0);
        _netArray[i]->setPartCount(1,0);
    }
    for(int i = 0; i < _cellNum; i++)
    {
        vector<int> list = _cellArray[i]->getNetList();
        if(!_cellArray[i]->getPart()){
            _partSize[0]++;
            for (int j = 0; j < list.size(); j++)
            {
                _netArray[list[j]]->incPartCount(0);
            }
        }
        else{
            _partSize[1]++;
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
    // cout<<_cutSize<<endl;
}
void Partitioner::add_cell(int i)
{
    if(_cellArray[i]->getLock()){
        return;
    }
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
            _bList[1][gain] = node;
        }
    }
}
void Partitioner::remove_cell(int i){
    if(_cellArray[i]->getLock()){
        return;
    }
    // cout<<"remove "<<i<<endl;
    Node *node = _cellArray[i]->getNode();
    if(!node->getPrev()){
        if(!node->getNext()){
            _bList[_cellArray[i]->getPart()].erase(_cellArray[i]->getGain());
        }else{
            _bList[_cellArray[i]->getPart()][_cellArray[i]->getGain()] = node->getNext();
            node->getNext()->setPrev(nullptr);
        }
    }else{
        if(node->getNext()){
            node->getPrev()->setNext(node->getNext());
            node->getNext()->setPrev(node->getPrev());
        }else{
            node->getPrev()->setNext(nullptr);
        }
    }
    node->setNext(nullptr);
    node->setPrev(nullptr);
}
void Partitioner::init()
{
    //init part size and count
    init_size_and_count();
    // init history
    _accGain = 0;
    _maxAccGain = 0;
    _moveNum = 0;
    _bestMoveNum = 0;
    _moveStack.clear();
    // init b_list
    _bList[0].clear();
    _bList[1].clear();
    //init gain and lock
    // go through every cell
    for (int i = 0; i < _cellNum; i++)
    {
        _cellArray[i]->setGain(0);
        _cellArray[i]->unlock();
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
        // std::cout << _cellArray[i]->getName() << " " << _cellArray[i]->getGain() << " " << _cellArray[i]->getPinNum() << " " << _cellArray[i]->getPart() << endl;
    }
}
void Partitioner::update_gain(int i)
{
    if(_cellArray[i]->getLock()){
        cout<<"Already locked, choose max is wrong"<<endl;
    }
    // cout<<"move "<<_cellArray[i]->getName()<<" gain "<<_cellArray[i]->getGain()<<endl;
    // update history
    _moveStack.push_back(i);
    _accGain += _cellArray[i]->getGain();
    _moveNum ++;
    if(_accGain > _maxAccGain){
        _maxAccGain = _accGain;
        _bestMoveNum = _moveNum;
    }
    // remove from bList and lock cell 
    remove_cell(i);
    _cellArray[i]->lock();
    //  update partSize
    if(!_cellArray[i]->getPart()){
        _partSize[0]-=1;
        _partSize[1]+=1;
    }else{
        _partSize[0]+=1;
        _partSize[1]-=1;
    }
     
    // check all nets
    vector<int> list = _cellArray[i]->getNetList();
    for (int j = 0; j < list.size(); j++)
    {
        // cout << _netArray[list[j]]->getName() << " " << _netArray[list[j]]->getPartCount(0) << " " << _netArray[list[j]]->getPartCount(1) << endl;
        vector<int> cells = _netArray[list[j]]->getCellList();

        // cell[i] in A(0)
        if (!_cellArray[i]->getPart())
        {
            // T(n)=0
            if (_netArray[list[j]]->getPartCount(1) == 0)
            {
                // cout << "T(n)=0" << endl;
                for (int k = 0; k < cells.size(); k++)
                {
                    remove_cell(cells[k]);
                    _cellArray[cells[k]]->incGain();
                    add_cell(cells[k]);
                    // cout << _cellArray[cells[k]]->getName() << "++" << endl;
                }
            }
            // T(n)=1
            else if (_netArray[list[j]]->getPartCount(1) == 1)
            {
                // cout << "T(n)=1" << endl;
                for (int k = 0; k < cells.size(); k++)
                {
                    if (_cellArray[cells[k]]->getPart())
                    {
                        remove_cell(cells[k]);
                        _cellArray[cells[k]]->decGain();
                        add_cell(cells[k]);
                        // cout << _cellArray[cells[k]]->getName() << "--" << endl;
                    }
                }
            }
            // lock cell[i]

            _netArray[list[j]]->incPartCount(1);
            _netArray[list[j]]->decPartCount(0);
            // F(n)=0

            if (_netArray[list[j]]->getPartCount(0) == 0)
            {
                // cout << "F(n)=0" << endl;
                for (int k = 0; k < cells.size(); k++)
                {
                    remove_cell(cells[k]);
                    _cellArray[cells[k]]->decGain();
                    add_cell(cells[k]);
                    // cout << _cellArray[cells[k]]->getName() << "--" << endl;
                }
            }
            // F(n)=1

            else if (_netArray[list[j]]->getPartCount(0) == 1)
            {
                // cout << "F(n)=1" << endl;
                for (int k = 0; k < cells.size(); k++)
                {
                    if (!_cellArray[cells[k]]->getPart())
                    {
                        remove_cell(cells[k]);
                        _cellArray[cells[k]]->incGain();
                        add_cell(cells[k]);
                        // cout << _cellArray[cells[k]]->getName() << "++" << endl;
                    }
                }
            }
        }
        // cell[i] in B(1)
        else
        {
            // T(n)=0
            if (_netArray[list[j]]->getPartCount(0) == 0)
            {
                // cout << "T(n)=0" << endl;
                for (int k = 0; k < cells.size(); k++)
                {
                    remove_cell(cells[k]);
                    _cellArray[cells[k]]->incGain();
                    add_cell(cells[k]);
                }
            }
            // T(n)=1
            else if (_netArray[list[j]]->getPartCount(0) == 1)
            {
                // cout << "T(n)=1" << endl;
                for (int k = 0; k < cells.size(); k++)
                {
                    if (!_cellArray[cells[k]]->getPart())
                    {
                        remove_cell(cells[k]);
                        _cellArray[cells[k]]->decGain();
                        add_cell(cells[k]);
                    }
                       
                }
            }
            // lock cell[i]
            _netArray[list[j]]->incPartCount(0);
            _netArray[list[j]]->decPartCount(1);
            // F(n)=0
            if (_netArray[list[j]]->getPartCount(1) == 0)
            {
                // cout << "F(n)=0" << endl;
                for (int k = 0; k < cells.size(); k++)
                {
                    remove_cell(cells[k]);
                    _cellArray[cells[k]]->decGain();
                    add_cell(cells[k]);
                }
            }
            // F(n)=1
            else if (_netArray[list[j]]->getPartCount(1) == 1)
            {
                // cout << "F(n)=1" << endl;
                for (int k = 0; k < cells.size(); k++)
                {
                    if (_cellArray[cells[k]]->getPart())
                    {
                        remove_cell(cells[k]);
                        _cellArray[cells[k]]->incGain();
                        add_cell(cells[k]);
                    }
                        
                }
            }
        }
    }
    // print gain
    // for (int i = 0; i < _cellArray.size(); i++)
    // {
    //     std::cout << _cellArray[i]->getName() << " " << _cellArray[i]->getGain() << endl;
    // }
}
void Partitioner::print_bList(int i){
    cout<<"bList "<<i<<endl;
    for(const auto& s : _bList[i]){
        
        cout<<s.first<<" ";
        Node *next = s.second;
        while (next)
        {
            cout<<_cellArray[next->getId()]->getName()<<" ";
            next = next->getNext();
        }
        
        cout<<endl;
    }
}
int Partitioner::choose_max(){
    map<int, Node *>::reverse_iterator itr0 = _bList[0].rbegin();
    map<int, Node *>::reverse_iterator itr1 = _bList[1].rbegin();
    if((_bList[0].empty())&&(_bList[1].empty())){
        // cout<<"Both empty"<<endl;
        return -1;
    }
    if(_bList[0].empty()){
        if(_partSize[1]-1 <min){
            // cout<<"A empty, B to small"<<endl;
            return -1;
        }
        // cout<<"A empty"<<endl;
        return (*itr1).second->getId();
    }
    if(_bList[1].empty()){
        if(_partSize[0]-1 <min){
            // cout<<"B empty, A to small"<<endl;
            return -1;
        }
        // cout<<"B empty"<<endl;
        return (*itr0).second->getId();
    }
    if(_partSize[0]-1 <min)
    {
        // cout<<"A to small"<<endl;
        return (*itr1).second->getId();
    }
        
    if(_partSize[1]-1 <min){
        // cout<<"B to small"<<endl;
        return (*itr0).second->getId();
    }
        
    int gain0 = _cellArray[(*itr0).second->getId()]->getGain();
    int gain1 = _cellArray[(*itr1).second->getId()]->getGain();
    if(gain0>gain1)
        return (*itr0).second->getId();
    if(gain0<gain1)
        return (*itr1).second->getId();
    if(_partSize[0]>=_partSize[1])
        return (*itr0).second->getId();
    return (*itr1).second->getId();
}
void Partitioner::partition()
{
    // initial partition
    init_part();
    init();

    // partition
    int iter = 0;
    do
    {
        init();
        while (1)
        {
            int index = choose_max();
            
            if(index < 0)
                break;
            // cout<<"move "<<_cellArray[index]->getName()<<" gain "<<_cellArray[index]->getGain()<<endl;
            update_gain(index);
        }
        for(int i=0;i<_bestMoveNum;i++){
            _cellArray[_moveStack[i]]->move();
        }
        cout<<"iter "<<iter<<" cutsize ";
        set_cutSize();
        iter++;
    } while (_maxAccGain>0);
    init();
    // cout<<"iter "<<iter<<" cutsize ";
    set_cutSize();
    init_size_and_count();
    // print_bList(0);
    // print_bList(1);
    
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
