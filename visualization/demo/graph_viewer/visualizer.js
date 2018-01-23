var nodes = {};

// Compute the distinct nodes from the links.
links.forEach(function(link) {
  link.source = nodes[link.source] || (nodes[link.source] = {name: link.source});
  link.target = nodes[link.target] || (nodes[link.target] = {name: link.target});
});

// var width = 960,
//     height = 500;

var force = d3.layout.force()
    .nodes(d3.values(nodes))
    .links(links)
    .size([width, height])
    .linkDistance(60)
    .charge(-200)
    .on("tick", tick)
    .start();

var svg = d3.select("#graph").append("svg")
    .attr("width", width)
    .attr("height", height);

// Per-type markers, as they don't inherit styles.
svg.append("defs").selectAll("marker")
    .data(["real"])
  .enter().append("marker")
    .attr("id", function(d) { return d; })
    .attr("viewBox", "0 -5 10 10")
    .attr("refX", 15)
    .attr("refY", 0)
    .attr("markerWidth", 6)
    .attr("markerHeight", 6)
    .attr("orient", "auto")
  .append("path")
    .attr("d", "M0,-5L10,0L0,5");

var path = svg.selectAll("text")
    .data(force.links())
  .enter().append("path")
    .attr("id", function(d) { return d.source.name + d.target.name; })
    .attr("class", function(d) {
      if (d.target.name.includes('$')) {
        return "link auxiliary";
      } else {
        return "link";
      }
    })
    .attr("marker-end", function(d) { return "url(#real)"; });

var edge_labels = svg.selectAll("text")
    .data(force.links())
  .enter().append("text")
    .style("font", "16px sans-serif");

edge_labels.append("textPath") //append a textPath to the text element
    .attr("xlink:href", function(d) {
      return "#" + d.source.name + d.target.name;
    }) //place the ID of the path here
    .style("text-anchor", "middle") //place the text halfway on the arc
    .attr("startOffset", "50%")
    .text(function(d) {
      return d.target.name[d.target.name.length - 1];
    });

var node = svg.selectAll(".node")
    .data(force.nodes())
  .enter().append("g")
    .attr("class", "node")
    .call(force.drag);

var circle = node.append("circle")
    .attr("r", 6)
    .style("fill", function(d) {
      if (d.name == "$".repeat(d.name.length)) {
        return "black";
      } else if (d.name.includes('$')) {
        if (d.name[0] != '$') {
          return "white";
        } else {
          return "gray";
        }
      } else {
        return "red";
      }
    })

var node_labels = node.append("text")
    .attr("x", 8)
    .attr("y", ".31em")
    .text(function(d) { return d.name; });

var margin = 10;
function tick() {
  node.each(function(d) {
    d.x = Math.max(margin, Math.min(width - margin, d.x));
    d.y = Math.max(margin, Math.min(height - margin, d.y));
  });
  node.attr("transform", transform);
  path.attr("d", linkArc);
}

function linkArc(d) {
  if (d.target.name == d.source.name) {
      // self-loop
      return "M" + (d.source.x - 0.01) + "," + d.source.y + "A" + 10 + "," + 12 + " 0 1,1 " + (d.target.x + 0.01) + "," + (d.target.y);
  }
  return "M" + d.source.x + "," + d.source.y + "L" + d.target.x + "," + d.target.y;
}

function transform(d) {
  return "translate(" + d.x + "," + d.y + ")";
}

