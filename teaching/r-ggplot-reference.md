# R & ggplot2 — Reference Sheet

A reference sheet I put together for myself: reusable `data.table` + `ggplot2` recipes for
figures that recur in applied empirical work. All data below is simulated and the numbers are
illustrative — each block is meant to be copied and adapted, not read as an analysis.

---

## 0. Setup — packages, palette, theme, saver

```r
library(data.table)   # dt[i, j, by]; ':=' adds/updates columns in place; .N = row count
library(ggplot2)
library(scales)       # axis labels: comma(), percent(), label_dollar(), label_number()

set.seed(42)

# one palette + one theme so every figure matches
pal <- c(blue = "#377EB8", red = "#E41A1C", green = "#4DAF4A", grey = "grey70")

theme_ref <- theme_classic(base_size = 13) +
  theme(axis.text     = element_text(colour = "black"),
        legend.position = "bottom",
        legend.title  = element_blank())

# consistent size; cairo_pdf embeds fonts cleanly
save_fig <- function(p, file, w = 8, h = 5)
  ggsave(file, p, width = w, height = h, device = cairo_pdf)
```

## A simulated panel used by the recipes below

```r
periods <- 1:60
event   <- 40L                                   # period the shock hits the treated unit
base    <- c(treated = 1000, ctrl_a = 800, ctrl_b = 1500, ctrl_c = 600)

dt <- CJ(group = names(base), period = periods)  # CJ = every group x period combination
dt[, value := base[group] *                      # index a named vector by the group column
             (1 - 0.002 * period) *              # mild common downward trend
             (1 + 0.03 * sin(2 * pi * period / 12)) *   # seasonality (12-period cycle)
             (1 + rnorm(.N, 0, 0.004))]          # noise
dt[group == "treated" & period >= event,         # treated: shock + gradual recovery
   value := value * (1 - 0.06 * exp(-(period - event) / 8))]
setkey(dt, group, period)
```

## 1. Fit a trend and extrapolate a counterfactual

Fit on the pre-event window only, then `predict` over *all* periods — the fit past the cutoff is
an out-of-sample projection. The gap between actual and projected is the object of interest.

```r
tr <- dt[group == "treated"]
tr[, season := factor(period %% 12)]             # seasonal dummies

# estimate on the pre-event window; predict everywhere = extrapolation
fit <- lm(log(value) ~ period + season, data = tr[period < event])
tr[, projected := exp(predict(fit, newdata = tr))]   # exp() undoes the log

ggplot(tr, aes(period, value)) +
  annotate("rect", xmin = event, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = "grey85", alpha = 0.5) +                       # shade the post-event region
  geom_ribbon(data = tr[period >= event],
              aes(ymin = value, ymax = projected),
              fill = pal[["red"]], alpha = 0.15) +               # the gap
  geom_line(colour = pal[["blue"]], linewidth = 1) +            # actual
  geom_line(aes(y = projected), colour = pal[["red"]], linetype = "22") +  # counterfactual
  labs(x = NULL, y = "Value") +
  theme_ref
```

## 2. Difference-in-differences, indexed

Index every series to 100 at a reference period so levels are comparable; the controls carry the
common trend, the treated series breaks from it.

```r
ref  <- event - 1L
dt[, idx := 100 * value / value[period == ref], by = group]      # rebase each series to 100
ctrl <- dt[group != "treated", .(idx = mean(idx)), by = period]  # control average

ggplot() +
  annotate("rect", xmin = event, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = "grey85", alpha = 0.5) +
  geom_hline(yintercept = 100, linetype = "dashed", colour = pal[["grey"]]) +
  geom_line(data = dt[group != "treated"], aes(period, idx, group = group),
            colour = pal[["grey"]], linewidth = 0.4) +           # individual controls, thin
  geom_line(data = ctrl, aes(period, idx),
            colour = pal[["green"]], linewidth = 1) +            # control average
  geom_line(data = dt[group == "treated"], aes(period, idx),
            colour = pal[["blue"]], linewidth = 1) +             # treated
  labs(x = NULL, y = "Index (= 100 at reference period)") +
  theme_ref
```

## 3. Event-study coefficients

Dynamic version of the DiD: leads and lags of the event, reference period `k = -1`. Flat before,
break at zero.

```r
es <- dt[period >= event - 12 & period <= event + 16]           # a window around the event
es[, k := fifelse(group == "treated", period - event, -1L)]     # controls pinned to reference
es[, k := factor(k, levels = c(-1L, setdiff(sort(unique(k)), -1L)))]  # -1 first = omitted

fit_es <- lm(log(value) ~ group + factor(period) + k, data = es) # unit FE + time FE + event time

co <- data.table(term = names(coef(fit_es)), b = coef(fit_es),
                 se = sqrt(diag(vcov(fit_es))))
co <- co[grepl("^k", term)][, kk := as.integer(sub("k", "", term))]
co <- rbind(co, data.table(term = "ref", b = 0, se = 0, kk = -1L))[order(kk)]  # add reference

ggplot(co, aes(kk, 100 * b)) +
  geom_hline(yintercept = 0, colour = pal[["grey"]]) +
  geom_vline(xintercept = -0.5, linetype = "dashed", colour = pal[["grey"]]) +
  geom_errorbar(aes(ymin = 100 * (b - 1.96 * se), ymax = 100 * (b + 1.96 * se)), width = 0) +
  geom_point(colour = pal[["blue"]], size = 2) +
  labs(x = "Periods since event", y = "Effect (%)") +
  theme_ref
```

## 4. Decompose a net change into two opposing flows

A signed bar chart: one flow above zero, the other below, so the net change and its two drivers
are visible at once.

```r
fl <- data.table(period = (event - 6L):(event + 10L))
fl[, inflow  := 40 + rnorm(.N, 0, 3)]
fl[, outflow := 45 + fifelse(period %in% event:(event + 2L), 60, 0) + rnorm(.N, 0, 3)]

fl <- melt(fl, id.vars = "period",                              # wide -> long for ggplot
           measure.vars = c("inflow", "outflow"),
           variable.name = "flow", value.name = "n")
fl[flow == "outflow", n := -n]                                  # draw outflow below zero

ggplot(fl, aes(period, n, fill = flow)) +
  annotate("rect", xmin = event - 0.5, xmax = event + 2.5, ymin = -Inf, ymax = Inf,
           fill = "grey85", alpha = 0.5) +
  geom_col(position = "identity", width = 0.8) +                # identity = overlay, not stack
  geom_hline(yintercept = 0, colour = pal[["grey"]]) +
  scale_fill_manual(values = c(inflow = pal[["green"]], outflow = pal[["red"]])) +
  labs(x = "Period", y = "Flow per period") +
  theme_ref
```

## 5. One-at-a-time sensitivity (tornado)

When a bottom-line number is a product of a few assumptions, vary each across a range holding the
others at their center. Bars sorted so the most consequential input is on top.

```r
bottom_line <- function(q = 50e3, margin = 12, life = 24) q * margin * life
central <- bottom_line()

sens <- data.table(
  input = c("Quantity", "Margin", "Lifetime"),
  lo = c(bottom_line(q = 40e3), bottom_line(margin = 10), bottom_line(life = 18)),
  hi = c(bottom_line(q = 60e3), bottom_line(margin = 14), bottom_line(life = 30)))
sens[, span := hi - lo]
setorder(sens, span)                                            # widest range last = plotted on top
sens[, input := factor(input, levels = input)]                  # lock the row order

ggplot(sens, aes(y = input)) +
  geom_vline(xintercept = central, linetype = "dashed", colour = pal[["grey"]]) +
  geom_segment(aes(x = lo, xend = hi, yend = input),
               colour = pal[["blue"]], linewidth = 6, alpha = 0.85, lineend = "round") +
  scale_x_continuous(labels = label_dollar(scale = 1e-6, suffix = "M")) +
  labs(x = "Bottom-line estimate", y = NULL) +
  theme_ref
```

## 6. Handy ggplot2 / scales patterns

```r
# --- axis labels (scales) ---
scale_y_continuous(labels = comma)                              # 1,000,000
scale_y_continuous(labels = label_number(scale = 1e-6, suffix = "M"))
scale_y_continuous(labels = percent)                           # 0.05 -> 5%
scale_x_continuous(labels = label_dollar())                    # $1,000

# --- shade a region / mark a point in time ---
annotate("rect", xmin = a, xmax = b, ymin = -Inf, ymax = Inf, fill = "grey85", alpha = 0.5)
geom_vline(xintercept = x0, linetype = "dashed")

# --- label series directly instead of a legend ---
annotate("text", x = x1, y = y1, label = "Treated", colour = pal[["blue"]], hjust = 0)

# --- reshape (data.table) ---
long <- melt(wide, id.vars = "id", measure.vars = c("a", "b"))  # wide -> long
wide <- dcast(long, id ~ variable)                              # long -> wide

# --- save ---
ggsave("fig.pdf", p, width = 8, height = 5, device = cairo_pdf)
```
