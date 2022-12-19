from django.contrib import admin
from django.urls import include, re_path
from django.conf import settings
from django.conf.urls.static import static
from normSeq import views

urlpatterns = [
    re_path('admin/', admin.site.urls),
    re_path(r'^$', views.index, name='home'),
    re_path(r'^query/', include ('query.urls')),
    re_path(r'^results/', include ('results.urls')),
    re_path(r'^status/', include ('status.urls')),
]+ static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)


