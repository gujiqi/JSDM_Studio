# UI startup fix

This build fixes the startup error:

```text
Error in b("Project") : could not find function "b"
```

Cause:
`b()` is not a Shiny/htmltools function. Bold text must use `tags$b()`.

Fix:
All bare `b("...")` calls were replaced with `tags$b("...")`, and a defensive alias `b <- tags$b` was added.
