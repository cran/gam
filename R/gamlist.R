"Gamlist" <-
function(...)
{
	gl <- list(...)
	oldClass(gl) <- c("Gamlist", "glmlist")
	gl
}
