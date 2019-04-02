function add_to_table(list, table_id) {
    if (list.length < 1)
        return;

    var columns = []

    for (var column_name in list[0]) {
        columns.push(column_name);
    }

    var tbl = document.getElementById(table_id);

    var tr = tbl.insertRow();

    tr.appendChild(document.createElement('th'))
      .appendChild(document.createTextNode("#"));

    for (var j = 0; j < columns.length; j++) {
        tr.appendChild(document.createElement('th'))
          .appendChild(document.createTextNode(columns[j]));
    }

    for (var i = 0; i < list.length; i++) {
        tr = tbl.insertRow();

        var dict = list[i];

        tr.appendChild(document.createElement('th'))
          .appendChild(document.createTextNode(i + 1));

        for (var j = 0; j < columns.length; j++) {
            tr.appendChild(document.createElement('td'))
              .appendChild(document.createTextNode(dict[columns[j]]));
        }
    }
}

add_to_table(samples, 'matched_samples');
