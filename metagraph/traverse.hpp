#ifndef __TRAVERSE_HPP__
#define __TRAVERSE_HPP__

#include <string>
#include <vector>

class DBG_succ;


namespace traverse {

    /**
     * Take the current graph content and return it in SQL
     * format (GA4GH Spec).
     *
     * We will perform one depth first search traversal of the graph. While we will record
     * one long reference string, we will output all sidepaths on the way.
     */
    void toSQL(DBG_succ *G, const std::vector<std::string> &fname,
                            const std::string &sqlfbase);

}

#endif // __TRAVERSE_HPP__
