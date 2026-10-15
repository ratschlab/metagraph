-- report_colour.lua: the colour marks of report.py (fmt_time) on the md -> pdf path (pandoc + xelatex).
-- report.md marks an execution time above report.TIME_MARK_YELLOW_MS / TIME_MARK_RED_MS as
-- <span class="t-yellow" style="background-color:#ffe866">N</span> (t-red, #ff8a80): pandoc reads the span as a
-- native Span, and for LaTeX output this filter turns the class into a \colorbox of the same tint (xcolor, which
-- pandoc's default LaTeX template loads) and defines the two colours in header-includes. Other output formats keep
-- the span (its inline style colours it where HTML renders). The text stays black: legible on both tints.
local colours = { ["t-yellow"] = { "suiteyellow", "FFE866" }, ["t-red"] = { "suitered", "FF8A80" } }

local function is_latex()
  return FORMAT:match("latex") ~= nil or FORMAT:match("beamer") ~= nil
end

function Span(el)
  if not is_latex() then
    return nil
  end
  for class, colour in pairs(colours) do
    if el.classes:includes(class) then
      -- \fboxsep 1pt (the default 3pt opens a visible gap before a comma after the box)
      local out = pandoc.Inlines({ pandoc.RawInline("latex", "{\\setlength{\\fboxsep}{1pt}\\colorbox{" .. colour[1] .. "}{") })
      for _, inline in ipairs(el.content) do
        out:insert(inline)
      end
      out:insert(pandoc.RawInline("latex", "}}"))
      return out
    end
  end
  return nil
end

function Meta(m)
  if not is_latex() then
    return nil
  end
  local defs = {}
  for _, colour in pairs(colours) do
    table.insert(defs, "\\definecolor{" .. colour[1] .. "}{HTML}{" .. colour[2] .. "}")
  end
  table.sort(defs)
  local block = pandoc.MetaBlocks({ pandoc.RawBlock("latex", table.concat(defs, "\n")) })
  local hi = m["header-includes"]
  if hi == nil then
    m["header-includes"] = pandoc.MetaList({ block })
  elseif hi.t == "MetaList" then
    hi:insert(block)
    m["header-includes"] = hi
  else
    m["header-includes"] = pandoc.MetaList({ hi, block })
  end
  return m
end
