
local action = _ACTION or ""
local examples = { "example1", "example2", "example3", "example4", "example5" }

solution "cpoly"
	location ( "build" )
	configurations { "Debug", "Release" }
	platforms {"native", "x64", "x32"}
  startproject "example4"
  
  for k,v in ipairs(examples) do
    project(tostring(v))
      kind "ConsoleApp"
      language "C++"
      files { "example/"..v..".c", "example/*.h", "src/*.h" }
      includedirs { "example", "src" }      
      targetdir "build"
          libdirs {"./src/GLFW/lib"}
     
      configuration { "linux" }
         links { "X11","Xrandr", "rt", "GL", "GLU", "pthread" }

      configuration { "windows" }
         links { "glu32","opengl32", "gdi32", "winmm", "user32" }

      configuration { "macosx" }
        links { "glfw3" }
        linkoptions { "-framework OpenGL", "-framework Cocoa", "-framework IOKit", "-framework CoreVideo" }

      configuration "Debug"
        defines { "DEBUG" }
        flags { "Symbols", "ExtraWarnings"}

      configuration "Release"
        defines { "NDEBUG" }
        flags { "Optimize", "ExtraWarnings"}    
  end