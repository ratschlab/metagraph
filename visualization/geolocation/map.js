/**
 * ---------------------------------------
 * This demo was created using amCharts 4.
 *
 * For more information visit:
 * https://www.amcharts.com/
 *
 * Documentation is available at:
 * https://www.amcharts.com/docs/v4/
 * ---------------------------------------
 */

// Themes begin
am4core.useTheme(am4themes_animated);
// Themes end

/**
 * Define SVG path for target icon
 */
var targetSVG = "M9,0C4.029,0,0,4.029,0,9s4.029,9,9,9s9-4.029,9-9S13.971,0,9,0z M9,15.93 c-3.83,0-6.93-3.1-6.93-6.93S5.17,2.07,9,2.07s6.93,3.1,6.93,6.93S12.83,15.93,9,15.93 M12.5,9c0,1.933-1.567,3.5-3.5,3.5S5.5,10.933,5.5,9S7.067,5.5,9,5.5 S12.5,7.067,12.5,9z";

// Create map instance
var chart = am4core.create("mapdiv", am4maps.MapChart);

// Set map definition
chart.geodata = am4geodata_worldLow;

// Set projection
// chart.projection = new am4maps.projections.Miller();
chart.projection = new am4maps.projections.Eckert6();

// Zoom control
chart.zoomControl = new am4maps.ZoomControl();

var homeButton = new am4core.Button();
homeButton.events.on("hit", function(){
  chart.goHome();
});

homeButton.icon = new am4core.Sprite();
homeButton.padding(7, 5, 7, 5);
homeButton.width = 20;
homeButton.icon.path = "M16,8 L14,8 L14,16 L10,16 L10,10 L6,10 L6,16 L2,16 L2,8 L0,8 L8,0 L16,8 Z M16,8";
homeButton.marginBottom = 10;
homeButton.parent = chart.zoomControl;
homeButton.insertBefore(chart.zoomControl.plusButton);

// Create map polygon series
var polygonSeries = chart.series.push(new am4maps.MapPolygonSeries());

/* Make map load polygon (like country names) data from GeoJSON */
polygonSeries.useGeodata = true;

/* Configure series */
var polygonTemplate = polygonSeries.mapPolygons.template;
polygonTemplate.applyOnClones = true;
polygonTemplate.togglable = true;
polygonTemplate.tooltipText = "{name}";
polygonTemplate.strokeWidth = 0.5;
polygonTemplate.strokeOpacity = 0.5;
// polygonTemplate.fill = chart.colors.getIndex(0).brighten(0.5).saturate(0.2);
polygonTemplate.fill = "#B4B4B7";
var lastSelected;
polygonTemplate.events.on("hit", function(ev) {
  if (lastSelected) {
    // This line serves multiple purposes:
    // 1. Clicking a country twice actually de-activates, the line below
    //    de-activates it in advance, so the toggle then re-activates, making it
    //    appear as if it was never de-activated to begin with.
    // 2. Previously activated countries should be de-activated.
    lastSelected.isActive = false;
  }
  ev.target.series.chart.zoomToMapObject(ev.target);
  if (lastSelected !== ev.target) {
    lastSelected = ev.target;
  }
})


/* Create selected and hover states and set alternative fill color */
var ss = polygonTemplate.states.create("active");
// ss.properties.fill = chart.colors.getIndex(5).brighten(0.0);
ss.properties.fill = chart.colors.getIndex(20).brighten(-0.3);

var hs = polygonTemplate.states.create("hover");
hs.properties.fill = chart.colors.getIndex(20).brighten(-0.1);

// Exclude Antartica
polygonSeries.exclude = ["AQ"];

// Small map
chart.smallMap = new am4maps.SmallMap();
// Re-position to top right (it defaults to bottom left)
chart.smallMap.align = "right";
chart.smallMap.valign = "top";
chart.smallMap.series.push(polygonSeries);

// create capital markers
var imageSeries = chart.series.push(new am4maps.MapImageSeries());

// define template
var imageSeriesTemplate = imageSeries.mapImages.template;
var circle = imageSeriesTemplate.createChild(am4core.Sprite);
circle.scale = 0.6;
circle.fill = chart.colors.getIndex(3).saturate(0.1).lighten(-0.5)
circle.path = targetSVG;
// what about scale...

// set propertyfields
imageSeriesTemplate.propertyFields.latitude = "latitude";
imageSeriesTemplate.propertyFields.longitude = "longitude";

imageSeriesTemplate.horizontalCenter = "middle";
imageSeriesTemplate.verticalCenter = "middle";
imageSeriesTemplate.align = "center";
imageSeriesTemplate.valign = "middle";
imageSeriesTemplate.width = 8;
imageSeriesTemplate.height = 8;
imageSeriesTemplate.nonScaling = true;
imageSeriesTemplate.tooltipText = "{sample}\n{city}";
imageSeriesTemplate.fill = am4core.color("#000");
imageSeriesTemplate.background.fillOpacity = 0;
imageSeriesTemplate.background.fill = "#fff";
imageSeriesTemplate.setStateOnChildren = true;
imageSeriesTemplate.states.create("hover");
imageSeries.data = [];
