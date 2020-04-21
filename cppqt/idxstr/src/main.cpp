#include <QCoreApplication>
#include <QFile>
#include <QTextStream>
#include <QDebug>



struct Link
{
    QString _name;
    qint64 _pos;
    bool operator<(const Link& other) { return _name < other._name; }
};



bool bubbleSort(QVector<Link>* list, const Link& link, int max)
{
    if (list->size()>=max)
    {
        return list->back()._name != link._name;
    }
    else
    {
        int i;
        for (i=0; i<list->size(); ++i)
        {
            if (link._name <= (*list)[i]._name)
            {
                list->insert(i,link);
                return true;
            }
        }
        list->insert(i,link);
        return true;
    }
}



int main(int argc, char *argv[])
{
    QCoreApplication a(argc,argv);
    QFile file(argv[1]);
    if (!file.open(QIODevice::ReadOnly|QIODevice::Text))
    {
        qDebug() << "Failed opening file.";
        return -1;
    }
    QTextStream in(&file);
    in.readLine();
    QVector<Link> lastList;
    QVector<Link> list;
    list.reserve(200000);
    while (!file.atEnd())
    {
        qDebug() << list.size();
        qint64 pos {in.pos()};
        QString line {in.readLine()};
        QStringList parts {line.split(" ")};
        Link link;
        link._name = parts[0];
        link._pos = pos;
        bubbleSort(&list,link,200000);
    }
    for (auto l: list)
    {
        qDebug() << l._name;
    }
    return 0;
}
